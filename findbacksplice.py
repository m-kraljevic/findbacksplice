#!/usr/bin/env python3
"""
FindBacksplice

Ran with 
    python blast_terminal_analysis.py -b <blast_xml_file> -i <probe_fasta> -g <genome_fasta>

Generates backsplice coordinates for use with circular RNA pipelines. Requires:
- the XML output of BLAST ran on a .fasta file of backsplice junction sequences/probes 
- the .fasta file of backsplice junction sequences/probes
- the .fasta file used to generate the BLAST database (i.e. a genome .fasta file)

Changing the default search length of 10,000 can be changed with the -s flag, for example a search of 100,000 bp can be accomplished with:
    python blast_terminal_analysis.py -b <blast_xml_file> -i <probe_fasta> -g <genome_fasta> -s 100000

"""

import sys
import argparse
import pandas as pd
from Bio.Blast import NCBIXML
from Bio import SearchIO, SeqIO

def parse_blast_xml(xml_file_path):
    # Parse BLAST XML and extract HSP data
    blast_data = {}
    
    with open(xml_file_path, 'r') as xml_file:
        blast_records = NCBIXML.parse(xml_file)
        
        for blast_record in blast_records:
            query_id = blast_record.query
            query_length = blast_record.query_length
            blast_data[query_id] = []
            
            for alignment in blast_record.alignments:
                hit_id = alignment.hit_id
                hit_def = alignment.hit_def
                
                for hsp in alignment.hsps:
                    hsp_data = {
                        'query_id': query_id,
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'query_length': query_length,
                        'hit_id': hit_id,
                        'hit_def': hit_def,
                        'hit_start': hsp.sbjct_start,
                        'hit_end': hsp.sbjct_end,
                        'strand': getattr(hsp, 'strand', ('Plus', 'Plus')),
                        'e_value': hsp.expect,
                        'bit_score': hsp.bits,
                        'identities': hsp.identities,
                        'align_length': hsp.align_length
                    }
                    blast_data[query_id].append(hsp_data)
    
    return blast_data

def blast_to_dataframe(blast_data):
    # Convert blast results to DataFrame
    all_hsps = []
    for query_id, hsps in blast_data.items():
        all_hsps.extend(hsps)
    return pd.DataFrame(all_hsps)

def filter_multi_hsp_hits(df):
    # Filter for hits with 2+ HSPs to be analyzed as potential pairs
    hit_hsp_counts = df.groupby(['query_id', 'hit_id']).size().reset_index(name='hsp_count')
    multi_hsp_hits = hit_hsp_counts[hit_hsp_counts['hsp_count'] >= 2]
    multi_hsp_pairs = set(zip(multi_hsp_hits['query_id'], multi_hsp_hits['hit_id']))
    mask = df.apply(lambda row: (row['query_id'], row['hit_id']) in multi_hsp_pairs, axis=1)
    return df[mask].copy()

def analyze_terminal_hsps(multi_hsp_df):
    # Analyze terminal (near start or near end) HSPs for hits with exactly 2 HSPs
    # this is to find potential backsplices based off only their BLAST output

    # Filter for hits with exactly 2 HSPs
    hit_counts = multi_hsp_df.groupby(['query_id', 'hit_id']).size().reset_index(name='hsp_count')
    two_hsp_hits = hit_counts[hit_counts['hsp_count'] == 2]
    two_hsp_pairs = set(zip(two_hsp_hits['query_id'], two_hsp_hits['hit_id']))
    two_hsp_mask = multi_hsp_df.apply(lambda row: (row['query_id'], row['hit_id']) in two_hsp_pairs, axis=1)
    two_hsp_df = multi_hsp_df[two_hsp_mask].copy()
    
    def format_strand(strand_tuple):
        """Convert strand tuple to simplified format"""
        if strand_tuple == ('Plus', 'Plus'):
            return '+'
        elif strand_tuple == ('Plus', 'Minus'):
            return '-'
        else:
            return str(strand_tuple)  # fallback for unexpected formats
    
    results = []
    
    for _, hit_info in two_hsp_hits.iterrows():
        query_id = hit_info['query_id']
        hit_id = hit_info['hit_id']
        
        hit_hsps = two_hsp_df[(two_hsp_df['query_id'] == query_id) & 
                             (two_hsp_df['hit_id'] == hit_id)].copy()
        
        if len(hit_hsps) != 2:
            continue
            
        query_length = hit_hsps.iloc[0]['query_length']
        hit_def = hit_hsps.iloc[0]['hit_def']
        
        # Find terminal HSPs
        c_terminal_hsp = hit_hsps[hit_hsps['query_end'] >= (query_length - 5)]
        n_terminal_hsp = hit_hsps[hit_hsps['query_start'] <= 5]
        
        # Final results to be reported by the tool
        result = {
            'circRNA': query_id,
            'chromosome': hit_def[0] if hit_def else '',
            'start': None,
            'end': None,
            'strand': None
        }
        
        if len(c_terminal_hsp) > 0:
            c_hsp = c_terminal_hsp.iloc[0]
            result['start'] = c_hsp['hit_start']
        
        if len(n_terminal_hsp) > 0:
            n_hsp = n_terminal_hsp.iloc[0]
            result['end'] = n_hsp['hit_end']
            result['strand'] = format_strand(n_hsp['strand'])
        
        results.append(result)
    
    return pd.DataFrame(results)

def load_sequences(fasta_file):
    # Load sequences from FASTA file into dict where the key 
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = record.seq
    return sequences

def backsplice_search(xml_file, probe_fasta, genome_fasta, processed_queries):
    # Perform a manual backsplice search using the genome file provided on all 
    # probes with a coordinate not generated from matching BLAST results

    # Load sequences
    circRNA_seqs = load_sequences(probe_fasta)
    genome = load_sequences(genome_fasta)
    blast_records = SearchIO.parse(xml_file, "blast-xml")
    
    results = []
    
    for blast_record in blast_records:
        for hit in blast_record:
            query_id = hit.hsps[0].query.id
            
            # Skip if already processed by pair matching
            if query_id in processed_queries:
                continue
            
            chromosome = hit.hsps[0].hit_id
            query_seq = circRNA_seqs.get(query_id)
            if not query_seq:
                continue
                
            found_hsp_pair = False
            
            # Collect HSP data
            query_starts = {}
            query_stops = {}
            hsp_objects = {}
            
            for i, hsp in enumerate(hit.hsps):
                query_start, query_end = hsp.query_range
                query_starts[query_start] = i
                query_stops[query_end] = i
                hsp_objects[i] = hsp
            
            query_length = len(query_seq)
            
            # Find HSPs near start and end
            near_start_hsps = [hsp_id for start, hsp_id in query_starts.items() if start < 10]
            near_end_hsps = [hsp_id for end, hsp_id in query_stops.items() if end > query_length - 10]
            
            # Search through near_start_hsps
            for start_hsp_id in near_start_hsps:
                if found_hsp_pair:
                    break
                    
                start_hsp = hsp_objects[start_hsp_id]
                start_hsp_range = start_hsp.query_range
                strand = start_hsp.hit_strand
                strand_symbol = "+" if strand == 1 else "-"
                
                search_seq = query_seq[start_hsp_range[1]:]
                
                if strand == -1:
                    search_seq = search_seq.reverse_complement()
                    search_start_coord = start_hsp.hit_range[1]
                    search_end_coord = start_hsp.hit_range[1] + 10000
                else:
                    search_start_coord = max(start_hsp.hit_range[0] - 10000, 0)
                    search_end_coord = start_hsp.hit_range[0]
                
                if search_seq and chromosome in genome:
                    search_chunk = genome[chromosome][search_start_coord:search_end_coord]
                    index = search_chunk.find(search_seq)
                    if index != -1:
                        found_hsp_pair = True
                        
                        if strand == -1:
                            start_coord = start_hsp.hit_range[0]
                            end_coord = search_start_coord + index + len(search_seq)
                        else:
                            start_coord = search_start_coord + index
                            end_coord = start_hsp.hit_range[1]
                        
                        result = {
                            'circRNA': query_id,
                            'chromosome': chromosome[0] if chromosome else '',
                            'start': start_coord,
                            'end': end_coord,
                            'strand': strand_symbol
                        }
                        results.append(result)
            
            # Search through near_end_hsps if not found
            if not found_hsp_pair:
                for end_hsp_id in near_end_hsps:
                    if found_hsp_pair:
                        break
                        
                    end_hsp = hsp_objects[end_hsp_id]
                    end_hsp_range = end_hsp.query_range
                    strand = end_hsp.hit_strand
                    strand_symbol = "+" if strand == 1 else "-"
                    
                    search_seq = query_seq[:end_hsp_range[0]]
                    
                    if strand == -1:
                        search_seq = search_seq.reverse_complement()
                        search_start_coord = max(end_hsp.hit_range[0] - 10000, 0)
                        search_end_coord = end_hsp.hit_range[0]
                    else:
                        search_start_coord = end_hsp.hit_range[1]
                        search_end_coord = end_hsp.hit_range[1] + 10000
                    
                    if search_seq and chromosome in genome:
                        search_chunk = genome[chromosome][search_start_coord:search_end_coord]
                        index = search_chunk.find(search_seq)
                        if index != -1:
                            found_hsp_pair = True
                            
                            if strand == -1:
                                start_coord = search_start_coord + index
                                end_coord = end_hsp.hit_range[1]
                            else:
                                start_coord = end_hsp.hit_range[0]
                                end_coord = search_start_coord + index
                            
                            result = {
                                'circRNA': query_id,
                                'chromosome': chromosome[0] if chromosome else '',
                                'start': start_coord,
                                'end': end_coord,
                                'strand': strand_symbol
                            }
                            results.append(result)
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description='BLAST Terminal HSP Analysis with Backsplice Detection')
    parser.add_argument('-b', '--blast', required=True, help='BLAST XML file')
    parser.add_argument('-i', '--input', required=True, help='Input probe FASTA file')
    parser.add_argument('-g', '--genome', required=True, help='Genome FASTA file')
    
    args = parser.parse_args()
    
    try:
        # First, run the backsplice generation via matching pairs of HSPs from BLAST
        blast_results = parse_blast_xml(args.blast)
        df = blast_to_dataframe(blast_results)
        multi_hsp_df = filter_multi_hsp_hits(df)
        terminal_results = analyze_terminal_hsps(multi_hsp_df)
        
        processed_queries = set(terminal_results['circRNA'].tolist())
        
        # Next, on all queries without a backsplice generated, perform a manual search
        backsplice_results = backsplice_search(args.blast, args.input, args.genome, processed_queries)
        
        # Finally, combine and save all results
        combined_results = pd.concat([terminal_results, backsplice_results], ignore_index=True)
        
        print(combined_results.to_csv(index=False), end='')
        
        print(f"# Terminal HSP results: {len(terminal_results)} queries", file=sys.stderr)
        print(f"# Backsplice search results: {len(backsplice_results)} queries", file=sys.stderr)
        print(f"# Total combined results: {len(combined_results)} queries", file=sys.stderr)
        print(f"# Output columns: circRNA, chromosome, start, end, strand", file=sys.stderr)
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
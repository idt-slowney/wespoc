#!/usr/bin/env python

from cgi import parse_header
import logging
from typing import Optional, List
from pathlib import Path
import re, os

import defopt
from collections import Counter, defaultdict
import multiprocessing as mp

import pandas as pd
import uuid
import statistics as stat
import subprocess
from subprocess import run


class Mpileup(object):
    """A mpileup class decribes a mpileup file

    Attributes:
        path    The path to the mpileup file

    Methods:
        fromBam        gen a mpileup file from bam file
        fromMpileup    load the mpileup file
        toTable        gen a parsed table from a mpileup
        toVcf          gen a vcf from a mpileup object by bcftools
        toErrorProfile gen an detailed error profile from a mpileup 
    """

    def __init__(self, path):
        assert os.path.exists(path), f"ERROR: {path} does not exist"
        self.path = path

    @classmethod
    def fromBam(cls, inbam:Path, maxDepth:int=10000, \
                samtools:Path='samtools', thread:int=4, \
                genome:Optional[Path]=None, inbed:Optional[Path]=None, \
                downSampleLevel:int=0, qualityThreshold:Optional[int]=None, \
                out:Optional[Path]=None):
        '''
        Read a bam file and return a mpileup object and file

        :param inbam:    the path to the inbam
            !!!inbam needs to be sorted and indexed
        :param maxDepth: the maxDepth when generating mpileup file.
            it is the -d option in samtools mpileup
        :param samtools: the path to the samtools
        :param thread:   the number of thread
        :param out:      the outpath of mpileup,
            When None, use the basename of the inbam file, the files would be
            saved in the same folder as the inbam file
        :param downSampleLevel: Down sample bam file before generating mpileup file
            When 0, there is no down sampling
        :param qualityThreshold: the quality threshold during mpileup
            When None, no minimal Q threshold is applied, Max: 40
        '''
        # check if in files exist
        for _ in [inbam, genome, inbed]:
            if _:
                assert os.path.exists(_), f"ERROR: {_} doesn't exist"
        # gen out mpileup file if out name not specified
        if not out:
            out = os.path.join(os.path.dirname(inbam), 
                               f'{Path(inbam).stem}.mpileup')

        mpileupCmd = f'{samtools} mpileup {inbam} -o {out}'
        # if we need to downsample the bam file
        if downSampleLevel:
            countOut = run(f'{samtools} view -c {inbam}'.split(' '), capture_output=True)
            totalReads = int(countOut.stdout.decode("utf-8").strip())
            dsFactor = downSampleLevel / totalReads
            # no ds needed
            if dsFactor > 1:
                logging.warning(f'DownSamplingLevel is higher than # of reads \
                                in {inbam}, No downsampling is done')
            else:
                mpileupCmd = f'{samtools} view -b -s {dsFactor} -t {thread} {inbam} | \
                            {samtools} mpileup - --output {out}'

        if maxDepth:
            mpileupCmd += f' -d {maxDepth}'
        if inbed:
            mpileupCmd += f' -l {inbed}'
        if genome:
            mpileupCmd += f' -f {genome}'
        if qualityThreshold:
            mpileupCmd += f' -Q {qualityThreshold}'
        # TODO: error out if mpileup fails
        subprocess.call(mpileupCmd, shell=True)

        return cls(out)

################################################################################
# this part turn mpile files into useful tables, vcfs, etc. for internal analyses
################################################################################    
    def toTable(self, output:Path, qualityCutOff:List[int]=[0, 20, 30], \
                strandness:bool=False, count:str='most_common', \
                ignore0covSites:bool=False, scoreOffSet:int=-33) -> None:
        ''' 
        Parses the mpileup and save a table describing base counts at each position
        Returns the path of the output table

        :param mpileup:       input mpileup file
        :param output:        path of the output table (a tab-delimited file)
        :param qualityCutOff: a list of quality cut offs
        :param strandness:    whether we separate the + and - strands counts
        :param count:         "all" or "most_common", 
            when all, output the count all the InDels
            when most_common, output the count the most common InDels
        :param ignore0covSites: ignore the sites with 0 coverage
        :param scoreOffSet:   the offset to calculate the Phred Q correctly
        '''
        logging.info(f"Parse the Mpileup file to a table: {self.path}")
        f = open(self.path, 'r')
        o = open(output, 'w')
        o.write(genHeader(qualityCutOff, strandness))

        for line in f:
            inline = mLine(line, strandness)
            # ignore the 0 coverage sites if specified by user
            if ignore0covSites:
                if inline.depth == 0:
                    continue
            # main logic
            # removing mapping info and indel, and record the InDels
            inline.removeMappingInfo()
            indelCounter = inline.parseOutIndels(countMode=count)
            # with differt Qual cut off, count the occurrence of ATGCN
            baseCounters = []
            for qual in qualityCutOff:
                basesAfterFilterQual = inline._filterBasesByQuality(qual, scoreOffSet)
                counts = inline.countBases(basesAfterFilterQual, strandness)
                baseCounters.append(BaseCounter(qual, counts, strandness))
            o.write(reportParsedMline(inline, baseCounters, indelCounter))
        f.close()
        o.close()

        logging.info("... Finished")
        return output
    
    def toTableParallel(self, output:Path, qualityCutOff:List[int]=[0, 20, 30], \
            strandness:bool=False, count:str='most_common', \
            ignore0covSites:bool=False, scoreOffSet:int=-33, thread:int=12, \
            tmpDir:str='/tmp') -> None:
        ''' 
        Parses the mpileup and save a table describing base counts at each position
        Returns the path of the output table

        :param mpileup:       input mpileup file
        :param output:        path of the output table (a tab-delimited file)
        :param qualityCutOff: a list of quality cut offs
        :param strandness:    whether we separate the + and - strands counts
        :param count:         "all" or "most_common", 
            when all, output the count all the InDels
            when most_common, output the count the most common InDels
        :param ignore0covSites: ignore the sites with 0 coverage
        :param scoreOffSet:   the offset to calculate the Phred Q correctly
        :param thread:        the number of thread to use
        '''
        logging.info(f"Parse the Mpileup file to a table: {self.path}")
        
        file_size = os.path.getsize(self.path)
        chunk_size = file_size // thread
        # print(f"{chunk_size=}") # python3.6 doesn't support this syntax
        
        # idea is to break the file into few chunks (# chunk = # threads)
        # find the start and end of each chunk and send 1 chunk to 1 thread
        # Arguments for each chunk 
        chunk_args = []
        with open(self.path, 'r') as f:
            def is_start_of_line(position):
                if position == 0:
                    return True
                # Check whether the previous character is EOL
                f.seek(position - 1)
                return f.read(1) == '\n'

            def get_next_line_position(position):
                # Read the current line till the end
                f.seek(position)
                f.readline()
                # Return a position after reading the line
                return f.tell()

            chunk_start = 0
            # Iterate over all chunks and construct arguments for `process_chunk`
            while chunk_start < file_size:
                chunk_end = min(file_size, chunk_start + chunk_size)

                # Make sure the chunk ends at the beginning of the next line
                while not is_start_of_line(chunk_end):
                    chunk_end -= 1

                # Handle the case when a line is too long to fit the chunk size
                if chunk_start == chunk_end:
                    chunk_end = get_next_line_position(chunk_end)

                # Save `process_chunk` arguments
                args = (self.path, chunk_start, chunk_end, ignore0covSites, qualityCutOff, strandness, scoreOffSet, count, tmpDir)
                chunk_args.append(args)

                # Move to the next chunk
                chunk_start = chunk_end

        with mp.Pool(thread) as p:
            # Run chunks in parallel
            chunk_results = p.starmap(Mpileup.processMpileupChunk, chunk_args)

        # writ ethe chunk output to the final output and remove the chunk output
        o = open(output, 'w')
        o.write(genHeader(qualityCutOff, strandness))
        
        for chunk_result in chunk_results:
            with open(chunk_result, 'r') as f:
                o.write(f.read())
        o.close()
        # delete all the tmp files
        for chunk in chunk_results:
            print(chunk)
            try:
                os.remove(chunk)
            except Exception as e:
                print(e)
        logging.info("... Finished")
        return output
    
    @staticmethod
    def processMpileupChunk(mpileup, start, end, ignore0covSites, qualityCutOff, strandness, scoreOffSet, count, tmpDir):
        '''
        returns the output file path of the processed chunk
        '''
        
        output = os.path.join(tmpDir, f"{uuid.uuid4()}.tmp")
        o = open(output, 'w')
        with open(mpileup, 'r') as f:
            # Moving stream position to `chunk_start`
            f.seek(start)
            # get the countDict for this chunk
            countDict = {}
            # Read and process lines until `chunk_end`
            for line in f:
                start += len(line)
                if start > end:
                    break
                inline = mLine(line, strandness)
                # ignore the 0 coverage sites if specified by user
                if ignore0covSites:
                    if inline.depth == 0:
                        continue
                # main logic
                # removing mapping info and indel, and record the InDels
                inline.removeMappingInfo()
                indelCounter = inline.parseOutIndels(countMode=count)
                # with differt Qual cut off, count the occurrence of ATGCN
                baseCounters = []
                for qual in qualityCutOff:
                    basesAfterFilterQual = inline._filterBasesByQuality(qual, scoreOffSet)
                    counts = inline.countBases(basesAfterFilterQual, strandness)
                    baseCounters.append(BaseCounter(qual, counts, strandness))
                o.write(reportParsedMline(inline, baseCounters, indelCounter))
        
        o.close()
        return output
                            
        
    def toVcf(self):
        pass #call bcltools to gen vcf files, do we need it?


    def toErrorProfileTable(self, output:Optional[Path]=None, countInDels:bool=False, \
                            count:str='most_common', ignore0covSites:bool=False, 
                            af:Optional[float]=None):
        '''
        Saves an error profile table that details the error rate and error weight
        of different types of errors
        Returns the error profile table path

        :param countInDels:  True or False, decides if we count InDels
        :param output:       the Path to the output table file
            When None, use the file name stem + error.profile.xlsx
        '''
        # turn mpileup to the parsed table first
        parsedTable = os.path.join(os.path.dirname(self.path), 
                                   f'{Path(self.path).stem}.mpileup.parsed.table')
        self.toTable(parsedTable, count=count, ignore0covSites=ignore0covSites)
        Mpileup.table2ErrorProfile(parsedTable, output, countInDels, af)

    def toErrorProfileTableParallel(self, output:Optional[Path]=None, countInDels:bool=False, \
                            count:str='most_common', ignore0covSites:bool=False, 
                            af:Optional[float]=None, threads:int=4, tmpDir:str="/tmp"):
        '''
        Saves an error profile table that details the error rate and error weight
        of different types of errors
        Returns the error profile table path

        :param countInDels:  True or False, decides if we count InDels
        :param output:       the Path to the output table file
            When None, use the file name stem + error.profile.xlsx
        '''
        # turn mpileup to the parsed table first
        parsedTable = os.path.join(os.path.dirname(self.path), 
                                   f'{Path(self.path).stem}.mpileup.parsed.table')
        self.toTableParallel(parsedTable, count=count, ignore0covSites=ignore0covSites, tmpDir=tmpDir, thread=threads)
        Mpileup.table2ErrorProfileParallel(parsedTable, output, countInDels, af, threads)

    @staticmethod
    def table2ErrorProfile(parsedTable:Path, output:Optional[Path]=None, \
                           countInDels:bool=False, af:Optional[float]=None):
        '''
        
        '''
        # df = pd.read_csv(parsedTable, sep='\t')
        
        errorDict = defaultdict(lambda: defaultdict(int))
        iterbaseString = 'A T G C'
        if countInDels:
            iterbaseString += ' insert_count deletion_count'
        with open(parsedTable, 'r') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                line = line.strip().split('\t')
                row = dict(zip(header, line))
                if af:
                    if float(row["AF"]) > af:
                        continue
                
                try:
                    refBase = row["ref"].upper()
                except:
                    print('error')
                    continue
                for base in iterbaseString.split():
                    errorDict[refBase][base] += int(row.get(base))
                    
            
        totalBase, errorBase = 0, 0
        
        # print(errorDict)

        outputDict = {}
        outputDict['Total Insertion count'] = 0
        outputDict['Total Deletion count'] = 0
        outputDict['Total base count'] = 0
        outputDict['Total error count'] = 0
        for ref in 'A T G C'.split(' '):
            for alt in 'A T G C'.split(' '):
                outputDict[f'Total {ref}>{alt} count'] = 0

        for ref, altDict in errorDict.items():
            for alt, count in altDict.items():
                if ref in 'A T G C'.split(' ') and alt in 'A T G C':
                    outputDict[f'Total {ref}>{alt} count'] = int(count)
                    if ref != alt:
                        outputDict['Total error count'] += int(count)
                if alt == 'insert_count':
                    outputDict['Total Insertion count'] += int(count)
                    outputDict['Total error count'] += int(count)
                if alt == 'deletion_count':
                    outputDict['Total Deletion count'] += int(count)
                    outputDict['Total error count'] += int(count)
                outputDict['Total base count'] += int(count)
        
        if outputDict['Total base count'] != 0:
            outputDict['TotalErrorRate(excludeN)'] = outputDict['Total error count'] / outputDict['Total base count']
            outputDict['TotalInDelRate(excludeN)'] = (outputDict['Total Insertion count'] + outputDict['Total Deletion count'])/ outputDict['Total base count']
        else:
            outputDict['TotalErrorRate(excludeN)'] = "NA"
            outputDict['TotalInDelRate(excludeN)'] = "NA"
        
        outputDict2 = {}
        for key, val in outputDict.items():
            if key != 'Total base count':
                tempKey = key.split('Total ')[-1].split(' count')[0]

                if outputDict['Total base count'] != 0:
                    rate = float(val / outputDict['Total base count'])
                else:
                    rate = "NA"

                if tempKey == 'Insertion' or tempKey == 'Deletion':
                    outputDict2[f'{tempKey} rate'] = rate
                else:
                    if tempKey.split('>')[0] != tempKey.split('>')[-1]:
                        outputDict2[f'{tempKey} substitution rate'] = rate


        outputDict = {**outputDict, **outputDict2}
        # print(outputDict)

        
        # output the table file
        errorProfile = os.path.join(os.path.dirname(parsedTable), 
                                    f'{Path(parsedTable).stem}.error.profile.tsv')
        
        # print(outputDict)
        if output:
            errorProfile = str(output)
            
        with open(output, 'w') as f:
            for key, value in outputDict.items():
                f.write(f'{key}\t{value}\n')
        # df = pd.DataFrame.from_dict(outputDict, orient='index')
        # # print(df)
        # df.to_csv(output, sep='\t', header=False)
        
        return errorProfile
    
    @staticmethod
    def table2ErrorProfileParallel(parsedTable:Path, output:Optional[Path]=None, \
                                   countInDels:bool=False, af:Optional[float]=None, \
                                   threads:int=4):
        '''
        
        '''
        
        file_size = os.path.getsize(parsedTable)
        chunk_size = file_size // threads
        # print(f"{chunk_size=}") # python3.6 doesn't support this syntax
        
        # idea is to break the file into few chunks (# chunk = # threads)
        # find the start and end of each chunk and send 1 chunk to 1 thread
        # Arguments for each chunk 
        chunk_args = []
        with open(parsedTable, 'r') as f:
            def is_start_of_line(position):
                if position == 0:
                    return True
                # Check whether the previous character is EOL
                f.seek(position - 1)
                return f.read(1) == '\n'

            def get_next_line_position(position):
                # Read the current line till the end
                f.seek(position)
                f.readline()
                # Return a position after reading the line
                return f.tell()

            chunk_start = 0
            # Iterate over all chunks and construct arguments for `process_chunk`
            while chunk_start < file_size:
                chunk_end = min(file_size, chunk_start + chunk_size)

                # Make sure the chunk ends at the beginning of the next line
                while not is_start_of_line(chunk_end):
                    chunk_end -= 1

                # Handle the case when a line is too long to fit the chunk size
                if chunk_start == chunk_end:
                    chunk_end = get_next_line_position(chunk_end)

                # Save `process_chunk` arguments
                args = (parsedTable, chunk_start, chunk_end, countInDels, af)
                chunk_args.append(args)

                # Move to the next chunk
                chunk_start = chunk_end

        with mp.Pool(threads) as p:
            # Run chunks in parallel
            chunk_results = p.starmap(Mpileup._processChunksErrorRate, chunk_args)

        errorDict = {}
        # print(f"{chunk_results=}")
        # Combine chunk results into `results`
        for chunk_result in chunk_results:
            # print(chunk_result)
            errorDict = combineCountDict(errorDict, chunk_result)        

        totalBase, errorBase = 0, 0
        
        # print(errorDict)

        outputDict = {}
        outputDict['Total Insertion count'] = 0
        outputDict['Total Deletion count'] = 0
        outputDict['Total base count'] = 0
        outputDict['Total error count'] = 0
        for ref in 'A T G C'.split(' '):
            for alt in 'A T G C'.split(' '):
                outputDict[f'Total {ref}>{alt} count'] = 0

        for ref, altDict in errorDict.items():
            for alt, count in altDict.items():
                if ref in 'A T G C'.split(' ') and alt in 'A T G C':
                    outputDict[f'Total {ref}>{alt} count'] = int(count)
                    if ref != alt:
                        outputDict['Total error count'] += int(count)
                if alt == 'insert_count':
                    outputDict['Total Insertion count'] += int(count)
                    outputDict['Total error count'] += int(count)
                if alt == 'deletion_count':
                    outputDict['Total Deletion count'] += int(count)
                    outputDict['Total error count'] += int(count)
                outputDict['Total base count'] += int(count)
                
        if outputDict['Total base count'] == 0:
            return
        
        if outputDict['Total base count'] != 0:
            outputDict['TotalErrorRate(excludeN)'] = outputDict['Total error count'] / outputDict['Total base count']
            outputDict['TotalInDelRate(excludeN)'] = (outputDict['Total Insertion count'] + outputDict['Total Deletion count'])/ outputDict['Total base count']
        else:
            outputDict['TotalErrorRate(excludeN)'] = "NA"
            outputDict['TotalInDelRate(excludeN)'] = "NA"
        
        outputDict2 = {}
        for key, val in outputDict.items():
            if key != 'Total base count':
                tempKey = key.split('Total ')[-1].split(' count')[0]

                if outputDict['Total base count'] != 0:
                    rate = float(val / outputDict['Total base count'])
                else:
                    rate = "NA"

                if tempKey == 'Insertion' or tempKey == 'Deletion':
                    outputDict2[f'{tempKey} rate'] = rate
                else:
                    if tempKey.split('>')[0] != tempKey.split('>')[-1]:
                        outputDict2[f'{tempKey} substitution rate'] = rate


        outputDict = {**outputDict, **outputDict2}
        # return outputDict

                            
        # output the table file
        errorProfile = os.path.join(os.path.dirname(parsedTable), 
                                    f'{Path(parsedTable).stem}.error.profile.tsv')
        
        # print(outputDict)
        if output:
            errorProfile = str(output)
            
        with open(errorProfile, 'w') as f:
            for key, value in outputDict.items():
                f.write(f'{key}\t{value}\n')
        # df = pd.DataFrame.from_dict(outputDict, orient='index')
        # # print(df)
        # df.to_csv(output, sep='\t', header=False)
        
        return errorProfile

    @staticmethod
    def calcErrorProfileVariance(parsedTable:Path, errorProfile:Path,
                                 allMpileup:Optional[str]=None,
                                 output:Optional[Path]=None,
                                 afThreshold:Optional[float]=0.005):
        '''
        Calculates variance (mean, stdev) of error rate AF.
        If a raw mpileup and UMI mpileups are provided, the mean and stdev
        will exclude sites where:
            1. Sites where ref matches the alt in the raw mpileup.
            2. Sites where ref != alt even after UMI consensus. In this case,
               the site may not be an error resulting from the wetlab protocol
               but instead may be a variant.
        '''
        rawMpileup = ''
        umiMpileupList = []
        print(allMpileup)
        if allMpileup != None:
            for x in allMpileup.split(','):
                basename = x.split('/')[-1].lower()
                if 'consensus' in basename:
                    umiMpileupList.append(x)
                elif 'umiduplicatedmarked' in basename:
                    rawMpileup = x
                elif 'duplicatemarked' in basename and rawMpileup == '':
                    rawMpileup = x

        ## Get all positions in rawMpileup where there is already 
        ## no error before UMI correction, to not count them as errors
        rawRefMatchingPositions = {}
        if rawMpileup != '':
            with open(rawMpileup, 'r') as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    if line.strip() == '':
                        continue
                    line = line.strip().split('\t')
                    row = dict(zip(header, line))
                    if (row['ref'].lower() == row['alt'].lower()) or (row['alt'].lower() == 'none'):
                        key = ':'.join([row['chrom'], row['position']])
                        rawRefMatchingPositions[key] = None

        ## Get all sites that were not corrected by UMI analysis
        ## This is because these may be true variants, possibly not errors
        uncorrectedSites = set()
        for mpileupTable in umiMpileupList:
            currentSites = []
            with open(mpileupTable, 'r') as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    if line.strip() == '':
                        continue
                    line = line.strip().split('\t')
                    row = dict(zip(header, line))
                    key = ':'.join([row['chrom'], row['position']])
                    if key in rawRefMatchingPositions:
                        # Don't count instances where there was no error in raw
                        continue
                    if float(row['AF']) > 0 and (row['ref'].lower() != row['alt'].lower()):
                        # Error found, add to ignore list
                        currentSites.append(key)
            if len(uncorrectedSites) == 0:
                uncorrectedSites = set(currentSites)
            else:
                # Only keep sites that were not corrected by either UMI analysis
                uncorrectedSites = uncorrectedSites.intersection(set(currentSites))

        # Convert to dict for quick lookup
        uncorrectedSites = {x: None for x in list(uncorrectedSites)}

        # Get error rate AF list to calculate variance
        errorAFList = []
        lowErrorAFList = []
        with open(parsedTable, 'r') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                if line.strip() == '':
                    continue
                line = line.strip().split('\t')
                row = dict(zip(header, line))
                key = ':'.join([row['chrom'], row['position']])
                af = float(row['AF'])

                if key in rawRefMatchingPositions:
                    errorAFList.append(0)
                elif key not in uncorrectedSites:
                    errorAFList.append(af)

        lowErrorAFList = [x for x in errorAFList if x <= afThreshold]

        outputDict = {}
        try:
            outputDict['TotalErrorRateMean(excludeN)'] = sum(errorAFList) / len(errorAFList)
        except:
            outputDict['TotalErrorRateMean(excludeN)'] = "NA"

        try:
            outputDict['TotalErrorRateStdDev(excludeN)'] = stat.stdev(errorAFList)
        except:
            outputDict['TotalErrorRateStdDev(excludeN)'] = "NA"

        try:
            outputDict['LowErrorRateStdDev(excludeN)'] = stat.stdev(lowErrorAFList)
        except:
            outputDict['LowErrorRateStdDev(excludeN)'] = "NA"
        
        with open(output, 'w') as o:
            with open(errorProfile, 'r') as i:
                for line in i:
                    o.write(line)
            for key, value in outputDict.items():
                o.write(f'{key}\t{value}\n')

    @staticmethod
    def _processChunksErrorRate(parsedTable:Path, chunk_start:int, chunk_end:int, countInDels, af):
        '''
        Go through a chunk of a file and count the numebr of errors
        '''
        errorDict = defaultdict(lambda: defaultdict(int))
        iterbaseString = 'A T G C'
        if countInDels:
            iterbaseString += ' insert_count deletion_count'
        with open(parsedTable, 'r') as f:
            header = f.readline().strip().split('\t')
            
            f.seek(chunk_start)
            for line in f:
                chunk_start += len(line)
                if chunk_start > chunk_end:
                    break

                line = line.strip().split('\t')
                
                if line == header:
                    # print("Skip")
                    continue
                    
                    
                row = dict(zip(header, line))
                
                if af:
                    if float(row["AF"]) >= af:
                        continue
                # print(row)
                
                try:
                    refBase = row["ref"].upper()
                except:
                    print('error')
                    continue
                for base in iterbaseString.split():
                    try:
                        errorDict[refBase][base] += int(row.get(base))
                    except Exception as e:
                        # print(e)
                        # an exception for the header line
                        pass
        
        # print(dict(errorDict))
        return dict(errorDict)

        
    @staticmethod
    def _groupErrorRates(errorRateDf):
        '''
        Returns a dataframe that group error rate of AT; GC together
        e.g. A>T and T>A is grouped to A>T|T>A
        ''' 
        # assumes the error rate df must have ATGC in index and columns
        result = {'A>T|T>A':0, 'A>G|T>C':0, 'A>C|T>G':0,
                  'G>A|C>T':0, 'G>T|C>A':0, 'G>C|C>G':0 }
        # TODO: clean up
        # awkard
        for refBase in 'ATGC':
            for altBase in 'ATGC':
                if refBase != altBase:
                    for key in result:
                        if f'{refBase}>{altBase}' in key:
                            rate = errorRateDf.loc[refBase, altBase]
                            result[key] += rate
        # TODO: when empty, what should we do?
        df = pd.DataFrame.from_dict(result, orient='index').rename(columns={0: 'rate'})

        return df#pd.DataFrame.from_dict(result).transpose()

################################################################################
# supporting classes and function
################################################################################
class mLine(object):
    def __init__(self, inputString, strandness):
        elements = inputString.strip().split('\t')
        try:
            self.chr = elements[0]
            self.pos = elements[1]
            self.ref = elements[2]
            self.depth = elements[3]
            self.bases = elements[4]
            self.quality = elements[5]
            if not strandness:
                self.ref = self.ref.upper()
        # cases there is no coverage, the bases and quality are empty
        except:
            self.chr = elements[0]
            self.pos = elements[1]
            self.ref = elements[2]
            self.depth = 0
            self.bases = ''
            self.quality = ''
            if not strandness:
                self.ref = self.ref.upper()

    def removeMappingInfo(self):
        ''' 
        Remove the ^. $ chars in the bases
        ^. represents the start of a read, the following char is the mapping Q
        $ represents the end of the quality
        '''
        newBases = re.sub(r'\^.', '', self.bases).replace('$','')
        self.bases = newBases

    def parseOutIndels(self, countMode='all'):
        ''' 
        Parse out the InDel part of the bases, input should have no mapping info
        remove the InDel patterns in the bases
        record the most frequent insert or deletion detected with its count
        Return: most common InDel and its count
        '''
        # init empty indel
        indelOut = InDelCounter(count=countMode)
        newBases = self.bases
        # deal with insertion
        if '+' in self.bases:
            # deal with all IUPAC nucleotide codes
            insertPattern = r"\+[0-9]+[ACGTNRYMKSWHBVDNacgtnrymkswhbvdn]+"
            insertions = re.findall(insertPattern, self.bases)
            cleanInserts = self._correctInDels(insertions)
            #remove inserts in bases
            for cInsert in set(cleanInserts.keys()):
                newBases = newBases.replace(cInsert, '')
            indelOut.updateInserts(cleanInserts, self.ref)
        # deal with deletion
        if '-' or '*' in self.bases:
            cleanDeletions = {}
            if '-' in self.bases:
                # print(self.bases)
                # deal with all IUPAC nucleotide codes
                delPattern = r"\-[0-9]+[ACGTNRYMKSWHBVDNacgtnrymkswhbvdn]+"
                deletions = re.findall(delPattern, self.bases)
                cleanDeletions = self._correctInDels(deletions)
                #remove inserts in bases
                # print(cleanDeletions)
                for cDel in set(cleanDeletions.keys()):
                    newBases = newBases.replace(cDel, '')
            # deal with *, weird case for *, it represent a deletion from previous position,
            # thus, it is represented in the previous position
            # we should ignore the *
            # weirdly, the quality of that base is in the quality strand
            if '*' in self.bases:
                pass
                # count = self.bases.count("*")
                # cleanDeletions['-0 '] = count
            # update the InDel counter
            if cleanDeletions != {}:
                indelOut.updateDeletions(cleanDeletions, self.ref)

        # assgin the new bases string without InDel info
        self.bases = newBases

        return indelOut

    def _correctInDels(self, inList):
        '''
        Correct the InDels based on the mpileup format
        '''
        if inList == []:
            return None
        
        result = defaultdict(int)
        for indel in inList:
            if indel not in result:
                indelLen = int(''.join(filter(str.isdigit, indel)))
                indelBases = indel.split(str(indelLen))[-1][:indelLen]
                modIndel = f"{indel[0]}{indelLen}{indelBases}"
                result[modIndel] += 1
            else:
                result[indel] += 1
        return result

    def _filterBasesByQuality(self, qualityCutOff, scoreOffSet):
        '''
        Filter out the Bases with lower Qualities if quality cut off exists
        '''
        # len of bases should be the same with len of quality at this step
        # print(self.bases)

        if len(self.bases) != len(self.quality):
            shorterLength = min(len(self.bases), len(self.quality))
            self.bases = self.bases[:shorterLength]
            self.quality = self.quality[:shorterLength]
            # f'Len of bases {len(self.bases)} and quality {len(self.quality)} are not consistent'
        # get the bases that are >= the quality threshold
        qualScores = [ord(i) + scoreOffSet for i in self.quality]
        filteredBases = ''.join([self.bases[_] \
            for _, qual in enumerate(qualScores) if qual >= qualityCutOff])
        return filteredBases


    def countBases(self, bases, strandness=False):
        '''
        Return a counter of differnt bases on a site based on the bases string
        '''
        if strandness:
            mixedBases = bases.replace('.', self.ref.upper()).replace(',', self.ref.lower())
            counter = Counter(mixedBases)
        else:
            upperBases = bases.replace('.', self.ref).replace(',', self.ref).upper()
            counter = Counter(upperBases)

        return counter
        
    def _getMaxOccuringChar(self, str, ASCII_SIZE = 256): 
        ''' 
        Returns the most frequent letters in a string, if input is empty, return a ""
        :param str: input string
        '''

        str = str.replace('.', '').replace(',', '')
        count = [0] * ASCII_SIZE 
        max = -1
        c = '' 

        for i in str: 
            count[ord(i)]+=1; 
        for i in str: 
            if max < count[ord(i)]: 
                max = count[ord(i)] 
                c = i 
        return c 

    def _mostFrequent(self, List):
        ''' returns the most frequent element in a list
        '''
        if List == []:
            return None
        return max(set(List), key=List.count)

def reportParsedMline(mline:mLine, baseCounters, \
                      indelCounter, sep:str='\t'):
    '''
    Returns a string in the output file after parsing the mpileup line

    :param mline: mline object
    :param baseCounter: baseCounter object about the bases
    :param indelCounter: indelCounter object bout the InDels
    :param sep: separator in output string
    '''

    out = []
    # first report position and depth
    out += [mline.chr, mline.pos, mline.ref, mline.depth]
    # second, report counter of bases at differnt Qual
    for baseCounter in baseCounters:
        counts = baseCounter.baseCounts
        if baseCounter.strandness:
            out += [counts.get(_) for _ in 'ATGCatgc']
        else:
            out += [counts.get(_) for _ in 'ATGC']

        AF = calcAf(mline, baseCounter, indelCounter)
        out.append(AF)
    # third, report indels
    out += [indelCounter.insert, indelCounter.inCount, \
               indelCounter.deletion, indelCounter.delCount]
    # add obAlt
    obAlt = getObAlt(mline.ref, baseCounters, indelCounter)
    out += [obAlt]
    return sep.join([str(_) for _ in out]) + '\n'

def calcAf(mLine, baseCounter, indelCounter):
    '''
    Calculate the max AF of alt allele
    '''
    counts = baseCounter.baseCounts
    ref = mLine.ref
    allCounts = sum([counts.get(_, 0) for _ in 'ATGCatgc']+ [indelCounter.inCount, indelCounter.delCount])
    ## TODO: this might be wrong
    if baseCounter.strandness:
        maxCount = max([counts.get(_, 0) for _ in 'ATGCatgc'.replace(ref.upper(), '').replace(ref.lower(), '')]+\
                [indelCounter.inCount, indelCounter.delCount])
    else:
        # when the ref is lower case or upper case
        maxCount = max([counts.get(_, 0) for _ in 'ATGCatgc'.replace(ref.upper(), '').replace(ref.lower(), '')]+\
                       [indelCounter.inCount, indelCounter.delCount])

    if allCounts != 0:
        AF = maxCount/ allCounts
    else:
        AF = 0
    return AF

def getObAlt(ref, baseCounters, indelCounter):
    '''
    Get the observed alt allele
    '''
    # use the counter without quality filtering
    counts = baseCounters[0].baseCounts
    tmp = 0
    alt = None
    for base in 'ATGC'.replace(ref, ''):
        if counts.get(base, 0) + counts.get(base.lower(), 0) > tmp:
            alt = base
            tmp = counts.get(base, 0) + counts.get(base.lower(), 0)
    if indelCounter.inCount > tmp:
        alt = indelCounter.insert
        tmp = indelCounter.inCount
    if indelCounter.delCount > tmp:
        alt = indelCounter.deletion
        tmp = indelCounter.delCount

    return alt

def genHeader(qualityCutOff, strandness, sep='\t'):
    '''
    Generate the first line for the output table
    '''
    out = []
    # first part
    out += "chrom position ref depth".split(' ')
    # bases part
    for qual in qualityCutOff:
        if strandness:
            for base in ['A', 'T', 'G', 'C', 'a', 't', 'g', 'c', 'AF']:
                if qual == 0:
                    out += [f'{base}']
                else:
                    out += [f'Q{qual}_{base}']
        else:
            for base in ['A', 'T', 'G', 'C', 'AF']:
                if qual == 0:
                    out += [f'{base}']
                else:
                    out += [f'Q{qual}_{base}']
    # indel part
    out += ['insert', 'insert_count', 'deletion', 'deletion_count']
    out += ['alt']

    return sep.join(out) + '\n'

def combineCountDict(a, b):
    '''
    combines 2 count dictionaries into 1
    '''
    result = {}
    # print(set(a) & set(b))
    for k in (set(a) | set(b)):
        # print(k)
        result[k] = combine_dicts(a.get(k, {}), b.get(k, {}))
        
    return result


def combine_dicts(a, b):
    '''
    adds 2 counter dictionaries together. e.g. {a: 1} + {a: 2} = {a: 3}
    '''
    return dict(Counter(a) + Counter(b))
class BaseCounter(object):

    def __init__(self, qualCutoff, counter, strandness):
        self.qualCutOff = qualCutoff
        self.strandness = strandness
        if strandness:
            self.baseCounts = {'A': counter['A'], 'T': counter['T'], \
                               'G': counter['G'], 'C': counter['C'], \
                               'a': counter['a'], 't': counter['t'], \
                               'g': counter['g'], 'c': counter['c']}
        else:
            self.baseCounts = {'A': counter['A'], 'T': counter['T'], \
                               'G': counter['G'], 'C': counter['C']}

class InDelCounter(object):
    def __init__(self, insert=None, insertCount=0, \
                 deletion=None, deletionCount=0, count='most_common'):
        self.insert = insert
        self.inCount = insertCount
        self.deletion = deletion
        self.delCount = deletionCount
        self.countMode = count

    def updateInserts(self, insertionDict, refBase):
        '''
        Update the output to the most common insert and its count 
        with the insertion dict
        '''
        tmpCounter = defaultdict(int)
        sum = 0
        for insert, count in insertionDict.items():
            insertLen = int(''.join(filter(str.isdigit, insert)))
            insertBase = refBase.upper() + insert.split(str(insertLen))[-1].upper()
            tmpCounter[insertBase] += count
            sum += count

        maxInsert, maxCount = max(tmpCounter.items(), key = lambda k : tmpCounter.get(k, -1))
        # maxInsert, maxCount = max(insertionDict.items(), key = lambda k : insertionDict.get(k, -1))
        # insertLen = int(''.join(filter(str.isdigit, maxInsert)))
        # outInsert = maxInsert.split(str(insertLen))[-1]
        self.insert = maxInsert
        self.inCount = maxCount
        if self.countMode == 'all':
            self.inCount = sum
        

    def updateDeletions(self, deletionDict, refBase):
        '''
        Update the output to the most common deletion and its count 
        with the insertion dict
        '''
        tmpCounter = defaultdict(int)
        sum = 0
        for delet, count in deletionDict.items():
            delLen = int(''.join(filter(str.isdigit, delet)))
            if delLen == 0:
                delBase = ''
            else:
                delBase = refBase.upper() + delet.split(str(delLen))[-1].upper()
            tmpCounter[delBase] += count
            sum += count

        maxDel, maxCount = max(tmpCounter.items(), key = lambda k : tmpCounter.get(k, -1))
        self.deletion = maxDel
        self.delCount = maxCount
        if self.countMode == 'all':
            self.delCount = sum

################################################################################
# cmd session of mpileup.py
################################################################################

def mpileupFromBam(*, inbam:Path, maxDepth:int=10000, \
                   samtools:Path='samtools', thread:int=4, \
                   genome:Optional[Path]=None, inbed:Optional[Path]=None, \
                   downSampleLevel:int=0, qualityThreshold:Optional[int]=None, \
                   out:Optional[Path]=None):
    '''Read a bam file and return a mpileup file, it can also downsample the bam file before generating the mpileup file

    :param inbam:    the path to the inbam
        !!!inbam needs to be sorted and indexed
    :param maxDepth: the maxDepth when generating mpileup file.
        it is the -d option in samtools mpileup
    :param samtools: the path to the samtools
    :param thread:   the number of thread used by samtools
    :param genome:   the genome to fetch the reference base during mpileup
    :param inbed:    the path of the inbed to limit the output position in mpileup
    :param out:      the outpath of mpileup,
        When None, use the basename of the inbam file, the files would be
        saved in the same folder as the inbam file
    :param downSampleLevel: # reads in the downsampled bam file
        When 0, there is no down sampling
    :param qualityThreshold: the quality threshold during mpileup
        When None, no minimal Q threshold is applied, Max: 40
    '''

    Mpileup.fromBam(inbam, maxDepth, samtools, thread, genome, 
                    inbed, downSampleLevel, qualityThreshold, out)


def mpileup2Table(*, mpileup:Path, output:Path, qualityCutOff:List[int]=[0, 20, 30], \
                strandness:bool=False, count:str='most_common', \
                ignore0covSites:bool=False, scoreOffSet:int=-33) -> None:
    '''
    Saves a the parsed table. For each position, we have depth, ref, alt, count of ATGC, inserts, and dels. Parses the mpileup and save a table describing base counts at each position

    :param mpileup:       input mpileup file
    :param output:        path of the output table (a tab-delimited file)
    :param qualityCutOff: a list of quality cut offs
    :param strandness:    whether we separate the + and - strands counts
    :param count:         "all" or "most_common", 
        when all, output the count of all the InDels
        when most_common, output the count of the most common InDels
    :param ignore0coverageSites: if a site has no coverage, ignore the site
    :param scoreOffSet:   the offset to calculate the Phred Q correctly
        -33 by default, we shouldn't change it
    '''
    mpileup = Mpileup(mpileup)
    mpileup.toTable(output, qualityCutOff, strandness, count, ignore0covSites, scoreOffSet)
    
    
def mpileup2TableParallel(*, mpileup:Path, output:Path, qualityCutOff:List[int]=[0, 20, 30], \
                strandness:bool=False, count:str='most_common', \
                ignore0covSites:bool=False, scoreOffSet:int=-33, \
                threads:int=12, tmpDir:str="/tmp") -> None:
    '''
    Saves a the parsed table. For each position, we have depth, ref, alt, count of ATGC, inserts, and dels. Parses the mpileup and save a table describing base counts at each position

    :param mpileup:       input mpileup file
    :param output:        path of the output table (a tab-delimited file)
    :param qualityCutOff: a list of quality cut offs
    :param strandness:    whether we separate the + and - strands counts
    :param count:         "all" or "most_common", 
        when all, output the count of all the InDels
        when most_common, output the count of the most common InDels
    :param ignore0coverageSites: if a site has no coverage, ignore the site
    :param scoreOffSet:   the offset to calculate the Phred Q correctly
        -33 by default, we shouldn't change it
    :param threads:       the number of threads to use
    '''
    mpileup = Mpileup(mpileup)
    mpileup.toTableParallel(output, qualityCutOff, strandness, count, ignore0covSites, scoreOffSet, threads, tmpDir)

def mpileup2ErrorProfile(*, mpileup:Path, output:Optional[Path]=None, \
                         countInDels:bool=False, count:str='all', \
                         ignore0covSites:bool=True, af:Optional[float]=None):
    '''
    Saves an error profile table that details the error rate and error weight of different types of errors.

    :param mpileup:         path of the input mpileup file
    :param count:           most_common or all, when all, AF equals (count of all alt allele) / depth, when most_common, AF is the AF of most common alt allele
    :param ignore0covSites: if we ignore the sites that has 0 coverage
    :param countInDels:     True or False, decides if we count InDels
    :param af:              a AF filter, when provided, the AF above this threshold would be ignored
    :param output:          the Path to the output table file
        When None, use the file name stem + error.profile.xlsx
    '''
    mpileup = Mpileup(mpileup)
    mpileup.toErrorProfileTable(output, count=count, ignore0covSites=ignore0covSites, \
                                countInDels=countInDels, af=af)

def calcErrorProfileVariance(*, table:Path, errorProfile:Path, allMpileup:Optional[str]=None, output:Optional[Path]=None, 
                            af:Optional[float]=0.005):
    '''
    Saves an error profile table that details the error rate and error weight of different types of errors.

    :param mpileup:      path of the input mpileup file
    :param errorProfile: path of error profile calculated by mpileup2ErrorProfile
    :param allMpileup:   comma-separated list of all mpileup files that will be used to exclude UMI consensus-corrected sites from variance calculation
    :param countInDels:  True or False, decides if we count InDels
    :param af:           a AF filter, when provided, the AF above this threshold would be ignored
    :param output:       the Path to the output table file
        When None, use the file name stem + error.profile.xlsx
    '''
    Mpileup.calcErrorProfileVariance(parsedTable=table, errorProfile=errorProfile, allMpileup=allMpileup, output=output, afThreshold=af)

def table2ErrorProfile(*, table:Path, output:Optional[Path]=None, countInDels:bool=True, \
                       af:Optional[float]=None):
    '''
    Saves an error profile table that details the error rate and error weight of different types of errors.

    :param mpileup:      path of the input mpileup file
    :param countInDels:  True or False, decides if we count InDels
    :param af:           a AF filter, when provided, the AF above this threshold would be ignored
    :param output:       the Path to the output table file
        When None, use the file name stem + error.profile.xlsx
    '''
    Mpileup.table2ErrorProfile(parsedTable=table, output=output, countInDels=countInDels, af=af)
    
def table2ErrorProfileParallel(*, table:Path, output:Optional[Path]=None, countInDels:bool=True, \
                       af:Optional[float]=None, threads:int=4):
    '''
    Saves an error profile table that details the error rate and error weight of different types of errors.

    :param mpileup:      path of the input mpileup file
    :param countInDels:  True or False, decides if we count InDels
    :param af:           a AF filter, when provided, the AF above this threshold would be ignored
    :param output:       the Path to the output table file
        When None, use the file name stem + error.profile.xlsx
    '''
    Mpileup.table2ErrorProfileParallel(parsedTable=table, output=output, countInDels=countInDels, af=af, threads=threads)


def bam2ErrorProfile(*, inbam:Path, inbed:Path, genome:Path,
                     maxDepth:int=10000, \
                     samtools:Path='samtools', thread:int=4, \
                     downSampleLevel:int=0, qualityThreshold:Optional[int]=None, \
                     ignore0covSites:bool=False, count:str='all', \
                     outMpileup:Optional[Path]=None, \
                     outErrorProfile:Optional[Path]=None, countInDels:bool=True, \
                     af:Optional[float]=None):
    '''
    Read a bam file and return a an error profile table, it can also downsample the bam file before generating the mpileup file

    :param inbam:       the path to the inbam
        !!!inbam needs to be sorted and indexed
    :param inbed:       the input bed file (usually a target bed file)
    :param maxDepth:    the maxDepth when generating mpileup file.
        it is the -d option in samtools mpileup
    :param samtools:    the path to the samtools
    :param thread:      the number of thread used by samtools
    :param outMpileup:  the outpath of mpileup,
        When None, use the basename of the inbam file, the files would be
        saved in the same folder as the inbam file
    :param downSampleLevel:     # reads in the downsampled bam file
        When 0, there is no down sampling
    :param qualityThreshold:    the quality threshold during mpileup
        When None, no minimal Q threshold is applied, Max: 40
    :param ignore0covSites:     if we ignore the sites that has 0 coverage
    :param count:       most_common or all, when all, AF equals (count of all alt allele) / depth, when most_common, AF is the AF of most common alt allele
    :param countInDels: True or False, decides if we count InDels
    :param af:          a AF filter, when provided, the AF above this threshold would be ignored
    :param output:      the Path to the output table file
        When None, use the file name stem + error.profile.xlsx
    '''

    mpileup = Mpileup.fromBam(inbam, maxDepth, samtools, thread, genome, 
                              inbed, downSampleLevel, qualityThreshold, outMpileup)
    mpileup.toErrorProfileTable(outErrorProfile, count=count, ignore0covSites=ignore0covSites, \
                                countInDels=countInDels, af=af)
    
def bam2ErrorProfileParallel(*, inbam:Path, inbed:Path, genome:Path,
                     maxDepth:int=10000, \
                     samtools:Path='samtools', thread:int=4, \
                     downSampleLevel:int=0, qualityThreshold:Optional[int]=None, \
                     ignore0covSites:bool=False, count:str='all', \
                     outMpileup:Optional[Path]=None, \
                     outErrorProfile:Optional[Path]=None, countInDels:bool=True, \
                     af:Optional[float]=None, tmpDir:str="/tmp"):
    '''
    Read a bam file and return a an error profile table, it can also downsample the bam file before generating the mpileup file

    :param inbam:       the path to the inbam
        !!!inbam needs to be sorted and indexed
    :param inbed:       the input bed file (usually a target bed file)
    :param maxDepth:    the maxDepth when generating mpileup file.
        it is the -d option in samtools mpileup
    :param samtools:    the path to the samtools
    :param thread:      the number of thread used by samtools
    :param outMpileup:  the outpath of mpileup,
        When None, use the basename of the inbam file, the files would be
        saved in the same folder as the inbam file
    :param downSampleLevel:     # reads in the downsampled bam file
        When 0, there is no down sampling
    :param qualityThreshold:    the quality threshold during mpileup
        When None, no minimal Q threshold is applied, Max: 40
    :param ignore0covSites:     if we ignore the sites that has 0 coverage
    :param count:       most_common or all, when all, AF equals (count of all alt allele) / depth, when most_common, AF is the AF of most common alt allele
    :param countInDels: True or False, decides if we count InDels
    :param af:          a AF filter, when provided, the AF above this threshold would be ignored
    :param output:      the Path to the output table file
        When None, use the file name stem + error.profile.xlsx
    '''

    mpileup = Mpileup.fromBam(inbam, maxDepth, samtools, thread, genome, 
                              inbed, downSampleLevel, qualityThreshold, outMpileup)
    mpileup.toErrorProfileTableParallel(outErrorProfile, count=count, ignore0covSites=ignore0covSites, \
                                countInDels=countInDels, af=af, threads=thread, tmpDir=tmpDir)


if __name__ == "__main__":
    defopt.run([mpileupFromBam, mpileup2Table, mpileup2ErrorProfile, \
                table2ErrorProfile, bam2ErrorProfile, mpileup2TableParallel, \
                table2ErrorProfileParallel, bam2ErrorProfileParallel,  calcErrorProfileVariance])
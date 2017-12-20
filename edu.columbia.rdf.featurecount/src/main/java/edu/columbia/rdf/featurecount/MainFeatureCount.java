package edu.columbia.rdf.featurecount;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.parsers.ParserConfigurationException;

import org.jebtk.bioinformatics.ext.samtools.SamUtils;
import org.jebtk.bioinformatics.gapsearch.GapSearch;
import org.jebtk.bioinformatics.gapsearch.SearchResults;
import org.jebtk.bioinformatics.genomic.Chromosome;
import org.jebtk.bioinformatics.genomic.GFF3Parser;
import org.jebtk.bioinformatics.genomic.GTB1Parser;
import org.jebtk.bioinformatics.genomic.GTB2Parser;
import org.jebtk.bioinformatics.genomic.GTBZParser;
import org.jebtk.bioinformatics.genomic.Gene;
import org.jebtk.bioinformatics.genomic.GeneParser;
import org.jebtk.bioinformatics.genomic.GeneType;
import org.jebtk.bioinformatics.genomic.Genes;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.Strand;
import org.jebtk.core.cli.CommandLineArg;
import org.jebtk.core.cli.CommandLineArgs;
import org.jebtk.core.cli.Options;
import org.jebtk.core.collections.CountMap;
import org.jebtk.core.collections.DefaultHashMap;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.DefaultTreeMapCreator;
import org.jebtk.core.collections.DoubleCountMap;
import org.jebtk.core.collections.HashSetCreator;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.IterTreeMap;
import org.jebtk.core.collections.TreeMapCreator;
import org.jebtk.core.collections.TreeSetCreator;
import org.jebtk.core.io.FileUtils;
import org.jebtk.core.io.PathUtils;
import org.jebtk.core.text.Join;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.xml.sax.SAXException;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamReader;


/**
 * Convert a sam file into a gene summary similar to htseq-count.
 * 
 * @author Antony Holmes Holmes
 *
 */
public class MainFeatureCount {
	private static final Logger LOG = 
			LoggerFactory.getLogger(MainFeatureCount.class);

	private static final int RECORD_N_NOTIFY = 100000;

	public static void main(String[] args) throws SAXException, IOException, ParserConfigurationException, ParseException {
		Options options = new Options();

		options.add('g', "gff", true);
		options.add('t', "gtb", true);
		options.add('e', "exclude-tag", true);
		options.add('l', "gene-level", true);
		options.add('o', "loc-frac");
		options.add('f', "gene-frac");
		options.add('s', "trans-frac");
		options.add('d', "detailed", true);
		options.add('m', "mapped");
		options.add('q', "quiet");
		options.add('u', "unmapped", true);
		options.add('b', "min-bp", true);

		boolean fracCountMode = false;

		boolean locFracMode = false;
		boolean geneFracMode = false;
		boolean transFracMode = false;
		boolean mapped = false;
		boolean quiet = false;

		// The fraction of cigar mappings that must contain the same gene for
		// a read to be kept
		double cigarKeepF = 0.6;

		Path gffFile = null;
		Path gtbFile = null;
		Path unmappedFile = null; //UNMAPPED_PATH;
		Path detailedFile = null; 
		int minBp = 10;

		Set<String> excludeTags = new HashSet<String>();
		Set<GeneType> levels = new HashSet<GeneType>();

		CommandLineArgs cmdArgs = CommandLineArgs.parse(options, args);

		GeneType level;

		
		for (CommandLineArg cmdArg : cmdArgs) {
			switch (cmdArg.getShortName()) {
			case 'g':
				gffFile = PathUtils.getPath(cmdArg.getValue());
				break;
			case 't':
				gtbFile = PathUtils.getPath(cmdArg.getValue());
				break;
			case 'e':
				excludeTags.add(cmdArg.getValue());
				break;
			case 'l':
				level = GeneType.parse(cmdArg.getValue());

				if (level != GeneType.UNDEFINED) {
					levels.add(level);
				}
				
				break;
			case 'b':
				minBp = cmdArg.getIntValue();
				break;
			case 'f':
				transFracMode = true;
				break;
			case 'o':
				locFracMode = true;
				break;
			case 's':
				geneFracMode = true;
				break;
			case 'd':
				detailedFile = PathUtils.getPath(cmdArg.getValue());
				break;
			case 'm':
				mapped = true;
				break;
			case 'q':
				quiet = true;
				break;
			case 'u':
				unmappedFile = PathUtils.getPath(cmdArg.getValue());
				break;
			case 'h':
			default:
				Options.printHelp(options);
				System.exit(0);
				break;
			}
		}

		if (levels.size() == 0) {
			levels.add(GeneType.EXON);
		}

		LOG.info("Excluding tags: {}...", excludeTags);


		//Map<String, Set<String>> transcriptSymbolMap = null;
		Map<String, Set<String>> symbolTranscriptMap = null;

		Genes genes = null;

		if (gffFile != null) {
			LOG.info("Loading gff from {}...", gffFile);

			GeneParser parser = new GFF3Parser()
					.setLevels(levels)
					.setKeepExons(false);

			genes = parser.parse(gffFile);
		}

		if (gtbFile != null) {
			LOG.info("Loading gtb from {}...", gtbFile);

			GeneParser parser;

			if (PathUtils.getName(gtbFile).contains("gtb2")) {
				parser = new GTB2Parser();
			} else if (PathUtils.getName(gtbFile).contains("gtbz")) {
				parser = new GTBZParser();
			} else {
				parser = new GTB1Parser();
			}

			parser = parser.setLevels(levels)
					.setKeepExons(false)
					.excludeByTag(excludeTags);

			genes = parser.parse(gtbFile);

			symbolTranscriptMap = 
					parser.idMap(gtbFile, "symbol", "transcript_id");
		}

		if (genes == null) {
			return;
		}

		fracCountMode = locFracMode || geneFracMode || transFracMode;

		//IterMap<String, Set<String>> readGeneMap =
		//		DefaultHashMap.create(new HashSetCreator<String>());

		//IterMap<String, Set<String>> geneTransMap =
		//		DefaultHashMap.create(new HashSetCreator<String>());

		IterMap<String, IterMap<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>>> readTransMap =
				DefaultTreeMap.create(
						DefaultTreeMapCreator.<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>>create(
								DefaultTreeMapCreator.<GenomicRegion, IterMap<String, Set<String>>>create(
										new DefaultTreeMapCreator<String, Set<String>>(new TreeSetCreator<String>()))));

		IterMap<String, IterMap<GenomicRegion, String>> cigarMap =
				DefaultTreeMap.create(new TreeMapCreator<GenomicRegion, String>());


		//IterMap<String, IterMap<GenomicRegion, Set<GenomicRegion>>> allCigarLocations =
		//		DefaultHashMap.create(
		//				DefaultHashMapCreator.<GenomicRegion, Set<GenomicRegion>>create(
		//						HashSetCreator.<GenomicRegion>create()));

		Set<GenomicRegion> allCigarLocations = new TreeSet<GenomicRegion>();

		//GapSearch<GFFGene> gapSearch = GFF.GFFToGapSearch(gff);

		//SAMRecordIterator iter = null;

		String bam = args[args.length - 1];

		Path bamFile = PathUtils.getPath(bam);

		//BufferedWriter countWriter = 
		//		FileUtils.newBufferedWriter(PathUtils.getPath(bam.replace(".bam", ".counts.txt")));

		if (!quiet) {
			LOG.info("Processing BAM from {}...", bamFile);

			int c = 1;

			//SAMTextHeaderCodec header = new SAMTextHeaderCodec();


			DoubleCountMap<String> tCountMap = DoubleCountMap.create();
			DoubleCountMap<String> gCountMap = DoubleCountMap.create();

			SamReader reader = SamUtils.newBamReader(bamFile);

			SAMFileWriter bamWriter = null;

			if (unmappedFile != null) {
				LOG.info("Writing unmapped reads to: {}", unmappedFile);

				bamWriter = SamUtils.newBamWriter(reader, unmappedFile);
			}

			try {
				for (final SAMRecord record : reader) {
					if (record.getReadUnmappedFlag()) {
						continue;
					}
					
					if (c % RECORD_N_NOTIFY == 0) {
						LOG.info("Processed {} records.", c);
					}

					++c;

					String rid = record.getReadName();

					String cigar = record.getCigarString();

					GenomicRegion region = SamUtils.getRegion(record);



					//int matches = record.getIntegerAttribute("NH");


					// Process what is left
					//GenomicRegion matchRegion = GenomicRegion.create(chr, start, end);

					allCigarLocations.clear();

					IterMap<GenomicRegion, SearchResults<Gene>> results = 
							processRegion(record, region, minBp, genes, allCigarLocations);

					if(!keep(record.getReadName(), region, results, allCigarLocations, cigarKeepF)) {
						if (bamWriter != null) {
							bamWriter.addAlignment(record);
						}

						continue;
					}




					//IterMap<String, Set<String>> tids =
					//		DefaultHashMap.create(new HashSetCreator<String>());

					cigarMap.get(rid).put(region, cigar);

					for (GenomicRegion cigarLocation : results) {
						SearchResults<Gene> search = results.get(cigarLocation);

						//System.err.println(rid + " " + results.size() + " " + search.size());

						for (GenomicRegion rr : search) {
							for (Gene gene : search.getValues(rr)) {
								String tid = gene.getTranscriptId();

								String symbol = gene.getSymbol(); //.toUpperCase();
								// for a given gene we know all of the transcripts
								// overlapping it
								//tids.get(symbol).add(tid);



								if (fracCountMode) {
									//transReadMap.get(tid).add(rid);
									readTransMap.get(rid).get(region).get(cigarLocation).get(symbol).add(tid);
									//readGeneMap.get(rid).add(symbol);

								} else {
									tCountMap.inc(tid);
									gCountMap.inc(symbol);
								}
							}
						}
					}

					
				}
			} finally {
				reader.close();

				if (bamWriter != null) {
					bamWriter.close();
				}
			}

			LOG.info("Processed {} records.", c);

			if (fracCountMode) {

				LOG.info("Counting...");

				for (String rid : readTransMap) {
					// How many genes this read is shared by

					IterMap<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>> locations = 
							readTransMap.get(rid);

					int locationCount;

					if (locFracMode) {
						locationCount = locations.size(); //readGeneMap.get(rid).size();
					} else {
						locationCount = 1;
					}

					// For each different location a read maps
					for (GenomicRegion location : locations) {

						// All the transcript ids of the elements it
						// overlaps
						IterMap<GenomicRegion, IterMap<String, Set<String>>> cigarLocations = 
								locations.get(location);

						//if (!keep2(rid, location, cigarLocations, allCigarLocations)) {
						//	continue;
						//}

						/*
						int geneCount;

						if (geneFracMode) {
							geneCount = countSymbols(cigarLocations);
						} else {
							geneCount = 1;
						}
						*/

						//CountMap<String> transCountsPerGene = null;

						//if(transFrac) {
						//	transCountsPerGene = countTranscriptsPerGene(cigarLocations);
						//}

						// How many cigar regions make up the read
						int cigarCount = cigarLocations.size();

						for (GenomicRegion cigarRegion : cigarLocations) {
							IterMap<String, Set<String>> symbols = 
									cigarLocations.get(cigarRegion);

							int geneCount;
							
							if (geneFracMode) {
								geneCount = symbols.size();
							} else {
								geneCount = 1;
							}
							
							for (String symbol : symbols) {
								Set<String> tids = symbols.get(symbol);

								int transCount;

								if (transFracMode) {
									// The read contributes 1 / transCount to the
									// count of this transcript since it overlaps
									// this many transcripts.
									transCount = tids.size(); //transCountsPerGene.get(symbol);
								} else {
									transCount = 1;
								}

								double f = f(locationCount, geneCount, cigarCount, transCount);

								for (String tid : tids) {
									//tCountMap.inc(tid, f);
									gCountMap.inc(symbol, f);
								}
							}
						}
					}
				}
			}


			output(genes, symbolTranscriptMap, gCountMap);
		}

		if (detailedFile != null) {
			detailedOutput(readTransMap, cigarMap, detailedFile);
		}

		if (mapped) {
			outputMapped(bam, bamFile, genes, minBp);
		}
	}

	private static void output(Genes genes, 
			Map<String, Set<String>> symbolTranscriptMap, 
			DoubleCountMap<String> gCountMap) {
		LOG.info("Writing...");

		Join.onTab("Gene", 
				"Chr", 
				"Start", 
				"End", 
				"Strand", 
				"Exon Width", 
				"Count", 
				"Transcript Ids")
		.println();


		Iterable<String> gids = genes.getSymbols();// CollectionUtils.sortKeys(gCountMap); //genes.getSymbols(); //

		for (String gid : gids) {
			Chromosome chr = null;
			Strand strand = null;
			//countWriter.write(name + "\t" + countMap.get(name));
			//countWriter.newLine();

			//String symbol = symbolTranscriptMap.get(gid).iterator().next();

			String transcripts = Join.onSemiColon()
					.values(symbolTranscriptMap.get(gid))
					.toString();

			Set<Integer> bases = new HashSet<Integer>();

			Iterable<Gene> exons = genes.getGenes(gid);

			int min = Integer.MAX_VALUE;
			int max = Integer.MIN_VALUE;

			// Count all the unique bases involved in each exon
			for (Gene exon : exons) {
				if (chr == null) {
					chr = exon.getChr();
				}

				if (strand == null) {
					strand = exon.getStrand();
				}

				min = Math.min(min, exon.getStart());
				max = Math.max(max, exon.getEnd());

				for (int i = exon.getStart(); i <= exon.getEnd(); ++i) {
					bases.add(i);
				}
			}

			Join.onTab(gid, 
					chr, 
					min, 
					max, 
					strand, 
					bases.size(), 
					gCountMap.get(gid), 
					transcripts)
			.println();
		}
	}


	private static boolean keep(String readId,
			GenomicRegion region,
			final IterMap<GenomicRegion, SearchResults<Gene>> results,
			final Set<GenomicRegion> allCigarLocations,
			double cigarKeepF) {


		int minC = (int)Math.ceil(allCigarLocations.size() * cigarKeepF);

		//System.err.println("cigars " + allCigarLocations + " " + minC);


		if (allCigarLocations.size() == 0) {
			return false;
		} else {
			// If there is more than one piece, check that they all share
			// a common gene

			Map<String, Set<GenomicRegion>> geneMap =
					DefaultHashMap.create(new HashSetCreator<GenomicRegion>());

			for (GenomicRegion cigarLocation : results) {
				SearchResults<Gene> se = results.get(cigarLocation);

				for (GenomicRegion re : se) {
					for (Gene gene : se.getValues(re)) {
						String symbol = gene.getSymbol();

						geneMap.get(symbol).add(cigarLocation);

						// If the gene is mapped in all of the cigar pieces,
						// we can assume that this is a spliced read
						if (geneMap.get(symbol).size() == minC) {
							return true;
						}
					}
				}
			}

			return false;
		}
	}

	private static final IterMap<GenomicRegion, SearchResults<Gene>> processRegion(SAMRecord record, 
			GenomicRegion region,
			int minBp,
			GapSearch<Gene> genes,
			Set<GenomicRegion> allCigarLocations) {
		//gapSearch.getOverlappingFeatures(region, minBp, results);

		// first lets look at the cigar

		Chromosome chr = region.getChr();
		int start = region.getStart();
		int end = start;
		Cigar cigar = record.getCigar();
		GenomicRegion r;

		IterMap<GenomicRegion, SearchResults<Gene>> ret =
				new IterTreeMap<GenomicRegion, SearchResults<Gene>>();

		for (CigarElement ce : cigar.getCigarElements()) {
			int l = ce.getLength();


			switch (ce.getOperator()) {
			case M:
				end = start + l - 1;

				r = GenomicRegion.create(chr, start, end);

				//System.err.println("CIGAR " + l + ce.getOperator() + " " + r);

				SearchResults<Gene> features = 
						genes.getOverlappingFeatures(r, minBp);

				// Only populate with interesting results
				if (features.size() > 0) {
					ret.put(r, features);
				}

				allCigarLocations.add(r);

				start = end + 1;

				break;
			case I:
				// Do nothing
				break;
			case D:
			case N:
				// Deletions mean we need to skip ahead in the reference
				start += l;
				break;
			default:
				break;

			}
		}

		return ret;
	}

	private static void detailedOutput(IterMap<String, IterMap<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>>> readTransMap,
			IterMap<String, IterMap<GenomicRegion, String>> samMap,
			Path file) throws IOException {
		LOG.info("Writing detailed output to {}...", file);

		BufferedWriter writer = FileUtils.newBufferedWriter(file);

		double sumF = 0;

		try {
			writer.write(Join.onTab().values("Read", "Location", "Width", "CIGAR Location", "CIGAR Width", "SAM", "Gene", "Transcript", "Location Count", "CIGAR Count", "Gene Count", "Transcript Count", "f").toString());
			writer.newLine();

			for (String rid : readTransMap) {
				// How many genes this read is shared by

				IterMap<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>> locations = 
						readTransMap.get(rid);

				int locationCount = locations.size(); //readGeneMap.get(rid).size();

				// For each different location a read maps
				for (GenomicRegion location : locations) {
					String sam = samMap.get(rid).get(location);

					// All the transcript ids of the elements it
					// overlaps
					IterMap<GenomicRegion, IterMap<String, Set<String>>> cigarLocations = 
							locations.get(location);

					//if (!keep2(rid, location, cigarLocations, allCigarLocations)) {
					//	continue;
					//}

					//String symbols = Join.onSemiColon().values(symbolMap).toString();

					int cigarCount = cigarLocations.size();
					
					//int geneCount = countSymbols(cigarLocations); //readGeneMap.get(rid).size();

					//CountMap<String> transCountsPerGene = 
					//		countTranscriptsPerGene(cigarLocations);

					

					for (GenomicRegion cigarLocation : cigarLocations) {
						IterMap<String, Set<String>> symbols = 
								cigarLocations.get(cigarLocation);

						int geneCount = symbols.size(); //countSymbols(cigarLocations);
						
						for (String symbol : symbols) {
							Set<String> tids = symbols.get(symbol);

							int transCount = tids.size();


							double f = f(locationCount, 
									cigarCount, 
									geneCount, 
									transCount);

							sumF += f;

							for (String tid : tids) {
								Join.onTab().values(rid, 
										location, 
										location.getLength(), 
										cigarLocation, 
										cigarLocation.getLength(), 
										sam, 
										symbol, 
										tid, 
										locationCount, 
										cigarCount, 
										geneCount, 
										transCount, 
										f)
								.println(writer);
							}
						}
					}
				}
			}
		} finally {
			writer.close();
		}

		LOG.info("sum(f) = {}", sumF);

		LOG.info("Finished.");
	}

	/**
	 * Calculate scaling factor for each read contribution.
	 * 
	 * @param locationCount
	 * @param geneCount
	 * @param cigarCount
	 * @param transCount
	 * @return
	 */
	private static double f(int locationCount, int cigarCount, int geneCount, int transCount) {
		int n = locationCount * cigarCount * geneCount * transCount;

		double f;

		if (n > 1) {
			f = 1.0 / n;
		} else {
			f = 1;
		}

		return f;
	}

	/**
	 * Count how many genes there are in total across each cigar region of
	 * a gene.
	 * 
	 * @param cigarLocations
	 * @return
	 */
	private static int countSymbols(IterMap<GenomicRegion, IterMap<String, Set<String>>> cigarLocations) {
		int ret = 0; //readGeneMap.get(rid).size();

		for (GenomicRegion cigarRegion : cigarLocations) {
			ret += cigarLocations.get(cigarRegion).size();
		}

		return ret;
	}

	private static CountMap<String> countTranscriptsPerGene(IterMap<GenomicRegion, IterMap<String, Set<String>>> cigarLocations) {
		CountMap<String> counts = new CountMap<String>();

		for (GenomicRegion cigarRegion : cigarLocations) {
			IterMap<String, Set<String>> symbols = cigarLocations.get(cigarRegion);

			for (String symbol : symbols) {
				Set<String> tids = symbols.get(symbol);

				counts.inc(symbol, tids.size());
			}
		}

		return counts;
	}

	private static void outputMapped(String bam, 
			Path bamFile,
			Genes genes,
			int minBp) throws IOException {
		
		String name = PathUtils.getName(PathUtils.getPwd());
		
		BufferedWriter transcriptomeUniqueWriter = 
				FileUtils.newBufferedWriter(PathUtils.getPath(name + ".ts-unique.sam"));

		BufferedWriter genomeUniqueWriter = 
				FileUtils.newBufferedWriter(PathUtils.getPath(name + ".ge-unique.sam"));

		BufferedWriter transcriptomeWriter = 
				FileUtils.newBufferedWriter(PathUtils.getPath(name + ".ts.sam"));

		BufferedWriter genomeWriter = 
				FileUtils.newBufferedWriter(PathUtils.getPath(name + ".ge.sam"));

		LOG.info("Processing BAM {} for mapping to {}...", bamFile, transcriptomeUniqueWriter);

		int c = 1;

		Set<GenomicRegion> allCigarLocations = new TreeSet<GenomicRegion>();

		SAMTextHeaderCodec header = new SAMTextHeaderCodec();

		try {
			SamReader reader = SamUtils.newBamReader(bamFile);

			header.encode(transcriptomeUniqueWriter, reader.getFileHeader());
			header.encode(genomeUniqueWriter, reader.getFileHeader());
			header.encode(transcriptomeWriter, reader.getFileHeader());
			header.encode(genomeWriter, reader.getFileHeader());

			for (final SAMRecord record : reader) {
				if (record.getReadUnmappedFlag()) {
					continue;
				}

				String sam = SamUtils.getSam(record);

				GenomicRegion region = SamUtils.getRegion(record);

				int matches = record.getIntegerAttribute("NH");

				allCigarLocations.clear();

				IterMap<GenomicRegion, SearchResults<Gene>> results =
						processRegion(record, region, minBp, genes, allCigarLocations);

				if (matches == 1) {
					if (results.size() > 0) {
						transcriptomeUniqueWriter.write(sam); // + "\t" + formatGenes(results));
						transcriptomeUniqueWriter.newLine();
					} else {
						//System.err.println("ge-unique " + sam);
						genomeUniqueWriter.write(sam);
						genomeUniqueWriter.newLine();
					}
				} else {
					if (results.size() > 0) {

						transcriptomeWriter.write(sam); // + "\t" + formatGenes(results));
						transcriptomeWriter.newLine();
					} else {
						//System.err.println("ge " + sam);
						genomeWriter.write(sam);
						genomeWriter.newLine();
					}
				}

				if (c % RECORD_N_NOTIFY == 0) {
					LOG.info("Processed {} records.", c);
				}

				++c;
			}

			reader.close();
		} finally {
			transcriptomeUniqueWriter.close();
			genomeUniqueWriter.close();
			transcriptomeWriter.close();
			genomeWriter.close();
		}
	}
}

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
import org.jebtk.bioinformatics.genomic.GTBParser;
import org.jebtk.bioinformatics.genomic.Gene;
import org.jebtk.bioinformatics.genomic.GeneParser;
import org.jebtk.bioinformatics.genomic.GeneType;
import org.jebtk.bioinformatics.genomic.Genes;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;


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
		options.add('b', "min-bp", true);

		boolean exonReads = false;

		boolean locFrac = false;
		boolean geneFrac = false;
		boolean transFrac = false;

		Path gffFile = null;
		Path gtbFile = null;

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
				transFrac = true;
				break;
			case 'o':
				locFrac = true;
				break;
			case 's':
				geneFrac = true;
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

			GeneParser parser = new GTBParser()
					.setLevels(levels)
					.setKeepExons(false)
					.excludeByTag(excludeTags);

			genes = parser.parse(gtbFile);

			symbolTranscriptMap = 
					parser.idMap(gtbFile, "symbol", "transcript_id");
		}

		if (genes == null) {
			return;
		}

		exonReads = locFrac || geneFrac || transFrac;

		//IterMap<String, Set<String>> readGeneMap =
		//		DefaultHashMap.create(new HashSetCreator<String>());

		//IterMap<String, Set<String>> geneTransMap =
		//		DefaultHashMap.create(new HashSetCreator<String>());

		IterMap<String, IterMap<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>>> readTransMap =
				DefaultTreeMap.create(
						DefaultTreeMapCreator.<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>>create(
								DefaultTreeMapCreator.<GenomicRegion, IterMap<String, Set<String>>>create(
										DefaultTreeMapCreator.<String, Set<String>>create(
												TreeSetCreator.<String>create()))));

		IterMap<String, IterMap<GenomicRegion, String>> cigarMap =
				DefaultTreeMap.create(TreeMapCreator.<GenomicRegion, String>create());


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

		LOG.info("Processing BAM from {}...", bamFile);

		int c = 1;

		//SAMTextHeaderCodec header = new SAMTextHeaderCodec();

		System.out.println(Join.onTab("Gene", "Transcript Ids", "Exon Width", "Count").toString());

		DoubleCountMap<String> tCountMap = DoubleCountMap.create();
		DoubleCountMap<String> gCountMap = DoubleCountMap.create();

		SamReader reader = 
				SamReaderFactory.makeDefault().open(bamFile.toFile());

		try {

			for (final SAMRecord record : reader) {
				if (record.getReadUnmappedFlag()) {
					continue;
				}

				String rid = record.getReadName();

				String cigar = record.getCigarString();

				GenomicRegion region = SamUtils.getRegion(record);



				int matches = record.getIntegerAttribute("NH");


				// Process what is left
				//GenomicRegion matchRegion = GenomicRegion.create(chr, start, end);

				allCigarLocations.clear();
				
				IterMap<GenomicRegion, SearchResults<Gene>> results = 
						processRegion(record, region, minBp, genes, allCigarLocations);

				if(!keep(record.getReadName(), region, results, allCigarLocations)) {
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

							String symbol = gene.getSymbol().toUpperCase();
							// for a given gene we know all of the transcripts
							// overlapping it
							//tids.get(symbol).add(tid);

							

							if (exonReads) {
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

				/*
				if (exonReads) {
					for (String tid1 : tids) {
						for (String tid2 : tids) {
							overlappingTransMap.get(tid1).add(tid2);
						}
					}
				}
				 */

				if (c % RECORD_N_NOTIFY == 0) {
					LOG.info("Processed {} records.", c);
				}

				++c;
			}

			if (exonReads) {
				
				LOG.info("Writing...", bamFile);
				
				for (String rid : readTransMap) {
					// How many genes this read is shared by

					IterMap<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>> locations = 
							readTransMap.get(rid);

					int locationCount;

					if (locFrac) {
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
						

						int geneCount;

						if (geneFrac) {
							geneCount = countSymbols(cigarLocations);
						} else {
							geneCount = 1;
						}

						//CountMap<String> transCountsPerGene = null;

						//if(transFrac) {
						//	transCountsPerGene = countTranscriptsPerGene(cigarLocations);
						//}

						// How many cigar regions make up the read
						int cigarCount = cigarLocations.size();

						for (GenomicRegion cigarRegion : cigarLocations) {
							IterMap<String, Set<String>> symbols = cigarLocations.get(cigarRegion);

							for (String symbol : symbols) {
								Set<String> tids = symbols.get(symbol);

								int transCount;

								if (transFrac) {
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


			Iterable<String> gids = genes.getSymbols();// CollectionUtils.sortKeys(gCountMap); //genes.getSymbols(); //

			for (String gid : gids) {
				//countWriter.write(name + "\t" + countMap.get(name));
				//countWriter.newLine();

				//String symbol = symbolTranscriptMap.get(gid).iterator().next();

				String transcripts = Join.onSemiColon().values(symbolTranscriptMap.get(gid)).toString();

				Set<Integer> bases = new HashSet<Integer>();

				Iterable<Gene> exons = genes.lookupBySymbol(gid);

				// Count all the unique bases involved in each exon
				for (Gene exon : exons) {

					GenomicRegion region = exon.getRegion();

					for (int i = region.getStart(); i <= region.getEnd(); ++i) {
						bases.add(i);
					}
				}


				System.out.println(Join.onTab(gid, transcripts, bases.size(), gCountMap.get(gid)).toString());
			}

		} finally {
			reader.close();
		}

		complexOutput(readTransMap, cigarMap, PathUtils.getPath("test2.txt"));
	}


	private static boolean keep(String readId,
			GenomicRegion region,
			final IterMap<GenomicRegion, SearchResults<Gene>> results,
			final Set<GenomicRegion> allCigarLocations) {
		
		System.err.println("cigars " + allCigarLocations);
		
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
						if (geneMap.get(symbol).size() == allCigarLocations.size()) {
							return true;
						}
					}
				}
			}

			return false;
		}
	}

	private static boolean keep2(String readId,
			GenomicRegion region,
			IterMap<GenomicRegion, IterMap<String, Set<String>>> cigarLocations,
			IterMap<String, IterMap<GenomicRegion, Set<GenomicRegion>>> allCigarLocations) {
		// All the locations of this read at this region
		Set<GenomicRegion> mappedCigarLocations = 
				allCigarLocations.get(readId).get(region);
		
		if (mappedCigarLocations.size() == 0) {
			return false;
		} else {
			// if there are multiple mapped locations, we must check that 
			// the same gene is found at each location

			// If there is more than one piece, check that they all share
			// a common gene

			Map<String, Set<GenomicRegion>> geneMap =
					DefaultHashMap.create(new HashSetCreator<GenomicRegion>());

			for (GenomicRegion cigarRegion : cigarLocations) {
				IterMap<String, Set<String>> symbols = cigarLocations.get(cigarRegion);

				for (String symbol : symbols) {
					//System.err.println("keep2 " + symbol + " " + cigarRegion + " " + mappedCigarLocations);
					
					geneMap.get(symbol).add(cigarRegion);

					// If the gene is mapped in all of the cigar pieces,
					// we can assume that this is a spliced read
					if (geneMap.get(symbol).size() == mappedCigarLocations.size()) { //cigarLocations.size()) {
						return true;
					}
				}
			}

			return false;
		}
	}

	private static final IterMap<GenomicRegion, SearchResults<Gene>> processRegion(SAMRecord record, 
			GenomicRegion region,
			int minBp,
			GapSearch<Gene> gapSearch,
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

				ret.put(r, gapSearch.getOverlappingFeatures(r, minBp));
				
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

	private static void complexOutput(IterMap<String, IterMap<GenomicRegion, IterMap<GenomicRegion, IterMap<String, Set<String>>>>> readTransMap,
			IterMap<String, IterMap<GenomicRegion, String>> samMap,
			Path file) throws IOException {
		LOG.info("Writing detailed output to {}...", file);

		BufferedWriter writer = FileUtils.newBufferedWriter(file);

		try {
			writer.write(Join.onTab().values("Read", "Location", "Width", "CIGAR Location", "CIGAR Width", "SAM", "Gene", "Transcript", "Locations", "Genes", "CIGARS", "Transcripts", "f").toString());
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

					int geneCount = countSymbols(cigarLocations); //readGeneMap.get(rid).size();

					//CountMap<String> transCountsPerGene = 
					//		countTranscriptsPerGene(cigarLocations);

					int cigarCount = cigarLocations.size();

					for (GenomicRegion cigarLocation : cigarLocations) {
						IterMap<String, Set<String>> symbols = cigarLocations.get(cigarLocation);

						for (String symbol : symbols) {
							Set<String> tids = symbols.get(symbol);

							int transCount = tids.size();


							double f = f(locationCount, geneCount, cigarCount, transCount);


							for (String tid : tids) {
								writer.write(Join.onTab().values(rid, location, location.getLength(), cigarLocation, cigarLocation.getLength(), sam, symbol, tid, locationCount, geneCount, cigarCount, transCount, f).toString());
								writer.newLine();
							}
						}
					}
				}
			}
		} finally {
			writer.close();
		}

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
	private static double f(int locationCount, int geneCount, int cigarCount, int transCount) {
		int n = locationCount * geneCount * cigarCount * transCount;

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
}
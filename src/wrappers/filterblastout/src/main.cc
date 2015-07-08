#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <map>

typedef struct Alignment {
	std::string qseqid = "", sseqid = "", sstrand = "", btop = "", qseq = "", sseq = "";
	int64_t qlen = 0, qstart = 0, qend = 0, slen = 0, sstart = 0, send = 0, score = 0, length = 0;
	float evalue = 0.0f, bitscore = 0.0f;
	std::string line;
} Alignment;

void filter_blast(std::string blast_out_path) {
	std::map<std::string, Alignment> top_alignments;

	std::ifstream infile(blast_out_path);
	std::string line;
	int64_t num_lines = 0;
	while (std::getline(infile, line)) {
		num_lines += 1;

		if ((num_lines % 1000) == 0) {
			fprintf (stderr, "\rProcessed %ld lines (kept: %ld)", num_lines, top_alignments.size());
		}

	    std::istringstream iss(line);

	   	// std::string qseqid, sseqid, sstrand, btop, qseq, sseq;
	   	// int64_t qlen = 0, qstart = 0, qend = 0, slen = 0, sstart = 0, send = 0, score = 0, length = 0;
	   	// float evalue = 0.0f, bitscore = 0.0f;
	   	Alignment alignment;
	   	alignment.line = line;
	    if (!(iss >> alignment.qseqid >> alignment.qlen >> alignment.qstart >> alignment.qend >> alignment.sseqid >> alignment.slen >> alignment.sstart >> alignment.send >> alignment.sstrand >> alignment.evalue >> alignment.bitscore >> alignment.score >> alignment.length >> alignment.btop >> alignment.qseq >> alignment.sseq))
	    	break;

	    // printf ("%s\t%ld\n", alignment.qseqid.c_str(), alignment.score);

		std::map<std::string, Alignment>::iterator it = top_alignments.find(alignment.qseqid);
		if (it != top_alignments.end()) {
			if (alignment.score > it->second.score) {
				it->second = alignment;
			}
		} else {
			top_alignments[alignment.qseqid] = alignment;
		}

		// if (num_lines > 10)
		// 	break;
	}

	fprintf (stderr, "\n");

	// printf ("\n");
	fprintf (stderr, "Number of unique qnames: %ld\n", top_alignments.size());

	int64_t num_registered_alignments = 0;
	for (std::map<std::string, Alignment>::iterator it = top_alignments.begin(); it != top_alignments.end(); it++) {
		num_registered_alignments += 1;
	    // printf ("[%ld] %s\t%ld\n", num_registered_alignments, it->second.qseqid.c_str(), it->second.score);
	    fprintf (stdout, "%s\n", it->second.line.c_str());
	}
}

void print_usage_and_exit(int argc, char **argv) {
	fprintf (stderr, "Tool for filtering alignments from BLAST's tabular format. The outfmt needs to be in specific order, hardcoded in the source of this program.\n");
	fprintf (stderr, "Outfmt: qseqid qlen qstart qend sseqid slen sstart send sstrand evalue bitscore score length btop qseq sseq\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage:\n");
	fprintf (stderr, "\t%s <alignments.out>\n", argv[0]);
	exit(0);
}

int main(int argc, char** argv) {
	if (argc != 2) {
		print_usage_and_exit(argc, argv);
	}

	std::string input_path = argv[1];

	// filter_blast("data/BLAST-v1.out");
	filter_blast(input_path);
	return 0;
}

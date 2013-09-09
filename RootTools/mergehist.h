#ifndef ROOT_MERGEHIST_H
#define ROOT_MERGEHIST_H

// author: Nils Krumnack (nils@krumnack.de), last change: 14 Oct 06

// This file includes a class that will merge the histograms from a
// series of file into one file, by adding them together.  It is
// intended to add together files from different runs from the same
// job, so it assumes all files have the same structure.  It does some
// sanity checks, but it is not totally failsafe.

// To use it:
//   MergeHist merger;
//   merger.input ("in1.root");
//   merger.input ("in2.root");
//   merger.input ("in3.root", 0.5);
//   merger.output ("out.root");
// This will add the histograms from in1.root, in2.root and in3.root, giving
// the last one half the weight of the others and write them out to out.root.


#include <string>
#include <vector>

class TDirectory;
class TString;

class MergeHist
{
  typedef std::vector<TDirectory*>::const_iterator iterator;

  // a list of all input files with associated weights.
private:
  std::vector<std::string> _input;
  std::vector<double> _weight;
public:
  inline void input (const std::string& input, double weight = 1) {
    _input.push_back (input); _weight.push_back (weight);};

  // write the combined histogams from the input files to an output file
public:
  void output (const std::string& output) const;

  // the actual functions for merging objects
private:
  void merge (TDirectory *out, const std::vector<TDirectory*>& dir) const;
  void merge_TH1 (TDirectory *out, const char *name,
		  const std::vector<TDirectory*>& dir) const;
  void merge_TDirectory (TDirectory *out, const char *name,
			 const std::vector<TDirectory*>& dir) const;
};

// this is the interface I spread in the past, it's use is depreceated now
void mergehist (const char *outName, int num, const TString *dirName);

#endif

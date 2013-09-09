//
// includes
//

#include "mergehist.h"

#include <iostream>
#include <TH1.h>
#include <TClass.h>
#include <TFile.h>

//
// methods of MergeHist
//

void MergeHist :: output (const std::string& output) const
{
  TFile out (output.c_str(), "RECREATE");
  TDirectory *old_dir = gDirectory;
  std::vector<TDirectory*> dir;

  try
  {
    for (std::vector<std::string>::const_iterator iter = _input.begin();
	 iter != _input.end(); ++ iter)
      dir.push_back (new TFile (iter->c_str()));

    merge (&out, dir);

    out.Write ();
  } catch (...)
  {
    for (unsigned i = 0; i < _input.size(); i ++)
      delete dir [i];
    old_dir->cd();
    throw;
  };
  for (unsigned i = 0; i < _input.size(); i ++)
    delete dir [i];
  old_dir->cd();
};

void MergeHist :: merge (TDirectory *out,
			 const std::vector<TDirectory*>& dir) const
{
  try
  {
    TObject *key;
    TIterator *iter;

    for (iter = dir [0]->GetListOfKeys ()->MakeIterator ();
	 (key = (*iter) ());)
    {
      const char *name = key->GetName ();

      if (!out->Get (name))
      {
	TObject *obj = dir [0]->Get (name);

	if (obj->InheritsFrom (TH1::Class ()))
	  merge_TH1 (out, name, dir);
	else if (obj->InheritsFrom (TDirectory::Class ()))
	  merge_TDirectory (out, name, dir);
	else
	  throw std::string (std::string ("Object type ") +
			     obj->IsA ()->GetName () + " unknown");
      };
    };
  } catch (std::string& s)
  {
    std::cerr << "Error: " << s << std::endl;
  };
};

void MergeHist :: merge_TH1 (TDirectory *out, const char *name,
			     const std::vector<TDirectory*>& dir) const
{
  std::vector<TH1*> hist;
  TH1 *result = NULL;

  try
  {
    for (iterator iter = dir.begin(); iter != dir.end(); ++ iter)
    {
      TObject *obj = (*iter)->Get (name);

      if (!obj) throw std::string ("object not found");
      if (!obj->InheritsFrom (TH1::Class ())) throw std::string ("object not a histogram");
      hist.push_back ((TH1 *) obj);
    };

    out->cd();
    result = (TH1 *) hist [0]->Clone ();
    result->Scale (_weight [0]);
    for (unsigned i = 1; i < dir.size(); i ++)
    {
      result->Add (hist [i], _weight [1]);
    };
  } catch (...)
  {
    delete result;
    throw;
  };
};

void MergeHist :: merge_TDirectory (TDirectory *out, const char *name,
				    const std::vector<TDirectory*>& dir) const
{
  std::vector<TDirectory*> dirs;
  TDirectory *result = NULL;

  for (unsigned i = 0; i < dir.size(); i ++)
  {
    TObject *obj = dir [i]->Get (name);

    if (!obj) throw std::string ("object not found");
    if (!obj->InheritsFrom (TDirectory::Class ())) throw std::string ("object not a directory");
    dirs.push_back ((TDirectory *) obj);
  };

  out->cd ();
  result = new TDirectory (dirs [0]->GetName (), dirs [0]->GetTitle ());
  merge (result, dirs);
};


void mergehist (const char *outName, int num, const TString *dirName)
{
  MergeHist merger;

  for (int i = 0; i < num; ++ i)
    merger.input (static_cast<const char*>(dirName[i]));
  merger.output (outName);
};

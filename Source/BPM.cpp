
#include <BPM.H>

#include <iostream>

using namespace std;

//
// Read springs from file.
//
// Each line represents a spring, format is: index, index, resting length, spring constant
//
void BPMMesh::InitFromAsciiFile(const string& file)
{
  // XXX: Use NReaders...

  ifstream ifs;
  ifs.open(file.c_str(), ios::in);

  if (!ifs.good())
    BoxLib::FileOpenFailed(file);

  int p1, p2;
  Real l, k;

  while (ifs >> p1 >> p2 >> l >> k) {
    Spring spring = { p1, p2, l, k };
    this->springs.push_back(spring);
    this->pmap.insert(pair<int,int>(p1, this->springs.size()-1));
    this->pmap.insert(pair<int,int>(p2, this->springs.size()-1));
  }

  cout << "BPMMesh::InitFromAsciiFile: " << this->springs.size() << " springs read." << endl;
}

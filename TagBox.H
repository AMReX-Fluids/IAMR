#ifndef _TagBox_H_
#define _TagBox_H_

//
// $Id: TagBox.H,v 1.5 1997-10-01 01:03:16 car Exp $
//

#include <IntVect.H>
#include <Box.H>
#include <Array.H>
#include <BaseFab.H>
#include <FabArray.H>
#include <BoxArray.H>
#include <BoxDomain.H>
#include <PArray.H>
#include <Geometry.H>

class TagBox
    : public BaseFab<int>
{
  public:
    enum TagVal{CLEAR=0, BUF, SET};

      // administrative
      // construct an invalid FAB with no memory
    TagBox();

    TagBox(const BOX &bx, int n = 1);

      // destructor
    virtual ~TagBox();
    
    void resize(const BOX& b) {BaseFab<int>::resize(b,1);}

    int& operator()(const IntVect& p)
    {
	return dptr[domain.index(p)];
    }
    const int& operator()(const IntVect& p) const
    {
	return dptr[domain.index(p)];
    }
 
      // construct and return a new tagbox in which the
      // coarsened cell is tagged of any of the corresponding
      // fine cells are tagged
    TagBox* coarsen(const IntVect & ratio);

      // mark neighbors of every tagged cell a distance nbuff away
      // only search interior for initial tagged points where nwid
      // is given as the width of the bndry region.
    void buffer(int nbuff, int nwid);

      // tag cells on intersect with src if corresponding src cell
      // is tagged.
    void merge(const TagBox& src);

      // add location of every tagged cell to IntVect array,
      // starting at given location
    int colate(Array<IntVect> &ar, int start) const;

#if 0
    void setVal(TagVal val) { BaseFab<int>::setVal(val,0); }
    void setVal(TagVal val, const BOX& bx) {
	BaseFab<int>::setVal(val,bx,0);
    }
#endif

    int numTags(const BOX& bx) const;
    int numTags() const;
  private:
    TagBox(const TagBox &);
    TagBox& operator = (const TagBox&);
};


class TagBoxArray : public FabArray<int, TagBox>
{
//protected:
    //PArray<TagBox> tb;
    int border;  // not to be confused with FabArray::n_grow which defines
		 // the invalid regions of the fabs in the FabArray
public:
    //TagBoxArray(PArrayPolicy _managed=PArrayManage)
	//: tb(_managed) {};
    //TagBoxArray(int len, PArrayPolicy _managed=PArrayManage)
	//: tb(len,_managed) {};
    //TagBoxArray(const BoxArray& bs, int _ngrow=0,
		//PArrayPolicy _managed = PArrayManage);

    TagBoxArray(const BoxArray& bs, int _ngrow=0);

    virtual ~TagBoxArray();

    int borderSize() { return border; }
    void buffer(int nbuf);
    void mergeUnique();

  // map tagged cells through a periodic boundary to other grids in TagBoxArray
  // cells which were outside domain are set to CLEAR
    void mapPeriodic( const Geometry& geom );

    void setVal(BoxDomain& bd, TagBox::TagVal val);
    void setVal(BoxArray& ba, TagBox::TagVal val);
    void coarsen(const IntVect & ratio);
    //int length() const { return tb.length(); }
    //TagBox& operator[](int i) {
	//return tb[i];
    //}
    //const TagBox& operator[](int i) const {
	//return tb[i];
    //}

    int numTags() const;

    Array<IntVect>* colate() const;
private:
    TagBoxArray();
    TagBoxArray(const TagBoxArray &);
    TagBoxArray& operator = (const TagBoxArray &);
};

#endif

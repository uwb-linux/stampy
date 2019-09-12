#include <iostream>

using std::min;
using std::max;

// Implements a traversal of the dynamic programming table along a diagonal band.
// Does not traverse position[0] == 0

class DiagonalBanding : Banding<2> {

 private:

  int iLen0, iLen1;
  int iWidth;
  Position pos;
  int iLastRow;

  int diagonal(int p1) { return ( 2*p1 - iLen1 + iLen0 + 1 ) / 2; }
    //return (p1 * iLen0 + iLen1/2) / iLen1;

 public:

  DiagonalBanding( int iLen0, int iLen1, int iWidth ) : iLen0(iLen0), iLen1(iLen1), iWidth(iWidth) {}

  Position& forwardIterator() {

    pos[0] = 1;
    pos[1] = 0;
    iLastRow = iLen0;
    return pos;

  }

  Position& backwardIterator() {

    pos[0] = iLen0;
    pos[1] = iLen1;
    iLastRow = 0;
    return pos;

  }

  bool hasNextForward() {

    if (pos[0] < iLen0 && pos[0] < (iWidth-1)/2 + max( diagonal(pos[1]+1)-1, diagonal(pos[1]))) {
      ++pos[0];
      return true;
    }
    if (pos[1]<iLen1) {
      ++pos[1];
      pos[0] = max(1, diagonal(pos[1]) - iWidth/2);
      return true;
    }
    return false;
  }

  bool hasNextBackward() {

    if (pos[0] > 1 && pos[0] > diagonal(pos[1]) - iWidth/2) {
      --pos[0];
      return true;
    }
    if (pos[1]>0) {
      --pos[1];
      pos[0] = max(1, min(iLen0, (iWidth-1)/2 + max( diagonal(pos[1]+1)-1, diagonal(pos[1]))));
      return true;
    }
    return false;
  }

  bool lastColumnEntry() {

    return pos[0] == iLastRow;

  }

  void warning() {

    cout << "Warning - out of bounds at position (" << pos[0] << "," << pos[1] << ")" << endl;

  }

};

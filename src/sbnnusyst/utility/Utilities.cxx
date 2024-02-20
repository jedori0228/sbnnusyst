#include "sbnnusyst/utility/Utilities.h"

namespace sbnnusyst{

int UniqueName(){
  static int N = 0;
  return N++;
}

} // END namespace sbnnusyst

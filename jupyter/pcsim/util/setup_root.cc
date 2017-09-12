// Simple ROOT script that imports the pcsim.util python module, which will
// automatically load the PCSIM libraries and necessary dependencies

#include <TPython.h>

void setup_root() { 
  TPython::Exec("import pcsim.util.setup_root"); 
  TPython::Exec("pcsim.util.setup_root.setup_root()");
}

# rucgraph

This is a library of cpp header files mainly used for graph computing tasks.


## How to use rucgraph

You can just download the whole repository into a rucgraph folder, and add the path of this folder when compiling cpp codes.

### An example of using rucgraph on a Linux server

1. Download the whole rucgraph folder, and save it on the server, e.g., the path of the saved folder is "root/rucgraph".

2. Given a cpp file named as "try.cpp", the contents of which are
```
using namespace std;

#include <data_structures/PairingHeapYS_with_offset.h>

int main()
{
	example_PairingHeapYS_with_offset();
}
```
, where "PairingHeapYS_with_offset.h" is a header file in rucgraph. In the terminal, compile the above file using the following command.
```
g++ -std=c++17 -I/root/rucgraph try.cpp a.out
./a.out
```



## Notes for updating files.

1. It is preferable to add an example in the end of each h file, for describing how to use the codes in the file.

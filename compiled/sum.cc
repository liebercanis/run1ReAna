#include "wsum.cc"
int main(int argc, char *argv[])
{
  cout << "executing " << argv[0] << endl;
  if (argc < 1)
  {
    printf(" usage: sum  <tag> <nevents> 0 is all  \n ");
    exit(0);
  }
  int nevents = 0;
  if (argc > 2)
    nevents = atoi(argv[2]);

  cout << argv[1] << endl;
  TString tag(argv[1]);

  new wsum(tag, nevents);
  exit(0);
}

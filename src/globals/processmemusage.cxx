#include<unistd.h>
/* #include<ios> */
/* #include<iostream> */
#include<fstream>
#include<string>
#include<limits>

//////////////////////////////////////////////////////////////////////////////
//
// processmemusage(double &, double &, double &) - takes three doubles by
// reference, attemps to read the system-dependent data for available and
// total memory as well as the system-dependent data for a process' resident
// set size, and return the results in KB.
//
// On failure, returns -1 on the failed value

extern"C" void processmemusage(double& memUsed, double& memAvail, double& memTotal)
{
   using std::ifstream;
   using std::string;

   memUsed  = -1;
   memAvail = -1;
   memTotal = -1;

   /* meminfo gives system totals */
   ifstream file("/proc/meminfo");

   file.ignore(18, ' ');
   file >> memTotal;
   file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

   // Skip 'MemFree:' line:
   file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

   file.ignore(18, ' ');
   file >> memAvail;
   file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

   file.close();

   /* /1* stat gives process totals *1/ */
   /* ifstream stat("/proc/self/stat"); */

   /* // dummy vars for leading entries in stat that we don't care about */
   /* // */
   /* string pid, comm, state, ppid, pgrp, session, tty_nr; */
   /* string tpgid, flags, minflt, cminflt, majflt, cmajflt; */
   /* string utime, stime, cutime, cstime, priority, nice; */
   /* string O, itrealvalue, starttime, vsize; */

   /* stat >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr */
   /*      >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt */
   /*      >> utime >> stime >> cutime >> cstime >> priority >> nice */
   /*      >> O >> itrealvalue >> starttime >> vsize >> memUsed; // don't care about the rest */

   /* stat.close(); */

   /* stat gives process totals */
   ifstream stat("/proc/self/statm");

   // dummy vars for leading entries in stat that we don't care about
   //
   string size;

   stat >> size >> memUsed; // don't care about the rest

   stat.close();

   double page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   memUsed = memUsed*page_size_kb;
}

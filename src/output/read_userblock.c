#include <string.h>
#include <stdio.h>
extern char userblock_start;
extern char userblock_end;
extern char userblock_size;

long get_userblock_size_(void) 
{
   //return (unsigned long)(&userblock_size);
   // Fixes mysterious bug occurring on some systems potentially due to the new GCC 7.3
   // where userblock_size is wrong though the symbol is correctly defined.
   // Since userblock_size = userblock_end - userblock_start , we just compute it on the fly.
   return (unsigned long)(&userblock_end-&userblock_start);
   
}

long get_inifile_size(char* filename) 
{
   FILE* fp = fopen(filename, "rb");
   fseek(fp,0,SEEK_END);
   long length=ftell(fp);
   fclose(fp);
   return length;
}

void insert_userblock(char* filename, char* inifilename, char* inifilename2) 
{
   FILE* fp = fopen(filename, "w");
   rewind(fp);
   fprintf(fp, "{[( START USERBLOCK )]}\n");

   // ini file
   fprintf(fp, "{[( INIFILE )]}\n");
   FILE* fini = fopen(inifilename, "rb");
   int c;
   do {
      c = fgetc (fini);
      if (c != EOF) fputc((char)c, fp);
   } while (c != EOF);
   fclose(fini);
   // DSMC ini file
   int lenIniFile2 = strlen(inifilename2);
   if (lenIniFile2>0) {
     fprintf(fp, "{[( DSMCFILE )]}\n");
     FILE* fini2 = fopen(inifilename2, "rb");
     int c;
     do {
        c = fgetc (fini2);
        if (c != EOF) fputc((char)c, fp);
     } while (c != EOF);
     fclose(fini2);
   }
   // compressed data ( as tar.xz )
   fprintf(fp, "{[( COMPRESSED )]}\n");
   fprintf(fp, "userblock.txt\n");      // filename
   fprintf(fp, "%ld\n", get_userblock_size_()); // filesize
   char* p = &userblock_start;
   while ( p != &userblock_end ) fputc(*p++, fp);
   fprintf(fp, "\n");

   fprintf(fp, "{[( END USERBLOCK )]}\n");
   fclose(fp);
}

void copy_userblock(char* outfilename, char* infilename)
{
   FILE* fout = fopen(outfilename,"rb+");
   FILE* fin  = fopen(infilename, "r");
   rewind(fout);
   int c;
   do {
      c = fgetc (fin);
      if (c != EOF) fputc((char) c, fout);
   } while (c != EOF);
   fclose(fin);
   fclose(fout);
}

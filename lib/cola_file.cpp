#include <cstdio>
#include <cassert>

#include "msg.h"
#include "cola_file.h"

using namespace std;

FILE* fp;
static int nhalo;

void cola_halo_file_open(const char filename[],
			 float* const boxsize)
{
  fp= fopen(filename, "r");
  if(fp == 0) {
    msg_printf(msg_fatal, "Error: unable to open cola fof file, %s\n",
	       filename);
    throw ColaFileError();
  }

  float params[3];
  int ret= fread(params, sizeof(float), 3, fp); assert(ret == 3);

  *boxsize= params[0];
  
  nhalo= 0;
}

int cola_halo_file_read_one(Halo* const h)
{
  int ret= fread(&h->nfof, sizeof(int), 1, fp); assert(ret == 1);
  if(h->nfof == 0) return 0;

  float f[3];
  
  ret= fread(h->x, sizeof(float), 3, fp); assert(ret == 3);
  ret= fread(h->v, sizeof(float), 3, fp); assert(ret == 3);
  ret= fread(f, sizeof(float), 3, fp);    assert(ret == 3);

  nhalo++;

  return 1;
}

void cola_halo_file_close()
{
  int nhalo_check= 0;
  int ret= fread(&nhalo_check, sizeof(int), 1, fp); assert(ret == 1);
  assert(nhalo == nhalo_check);

  ret= fclose(fp); assert(ret == 0);
}

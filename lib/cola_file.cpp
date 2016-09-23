#include <cstdio>
#include <cassert>

#include "msg.h"
#include "cola_file.h"

using namespace std;

FILE* fp_fof;
static int nhalo;

FILE* fp_part;
static int npart, npart_read;

void cola_halo_file_open(const char filename[],
			 float* const boxsize)
{
  fp_fof= fopen(filename, "r");
  if(fp_fof == 0) {
    msg_printf(msg_fatal, "Error: unable to open cola fof file, %s\n",
	       filename);
    throw ColaFileError();
  }

  float params[3];
  int ret= fread(params, sizeof(float), 3, fp_fof); assert(ret == 3);

  *boxsize= params[0];
  
  nhalo= 0;
}

int cola_halo_file_read_one(Halo* const h)
{
  int ret= fread(&h->nfof, sizeof(int), 1, fp_fof); assert(ret == 1);
  if(h->nfof == 0) return 0;

  float f[3];
  
  ret= fread(h->x, sizeof(float), 3, fp_fof); assert(ret == 3);
  ret= fread(h->v, sizeof(float), 3, fp_fof); assert(ret == 3);
  ret= fread(f, sizeof(float), 3, fp_fof);    assert(ret == 3);

  nhalo++;

  return 1;
}

void cola_halo_file_close()
{
  int nhalo_check= 0;
  int ret= fread(&nhalo_check, sizeof(int), 1, fp_fof); assert(ret == 1);
  assert(nhalo == nhalo_check);

  ret= fclose(fp_fof); assert(ret == 0);
}


void cola_part_file_open(const char filename[],
			 float* const boxsize, int* np)
{
  fp_part= fopen(filename, "r");
  if(fp_part == 0) {
    msg_printf(msg_fatal, "Error: unable to open cola particle file, %s\n",
	       filename);
    throw ColaFileError();
  }

  float header[6];
  int ret= fread(header, sizeof(float), 6, fp_part); assert(ret == 6);
  ret = fread(&npart, sizeof(int), 1, fp_part);
  npart_read= 0;

  *boxsize= header[0];
  *np= npart;
}

void cola_part_file_close()
{
  float buf[9];
  for(int i=npart_read; i<npart; ++i) {
    int ret= fread(buf, sizeof(float), 9, fp_part);
    assert(ret == 9);
  }

  int np_check;
  int ret= fread(&np_check, sizeof(int), 1, fp_part); assert(ret == 1);
  assert(np_check == npart);

  ret= fclose(fp_part); assert(ret == 0);
  npart= 0;
}


int cola_part_file_read_one(Halo* const h)
{
  if(npart == npart_read) return 0;

  float f[3];
  int ret= fread(h->x, sizeof(float), 3, fp_part); assert(ret == 3);
  ret= fread(h->v, sizeof(float), 3, fp_part); assert(ret == 3);
  ret= fread(f, sizeof(float), 3, fp_part); assert(ret == 3);

  npart_read++;

  return 1;
}

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

double dist(double a1[], double a2[]) {
  double dd;
  dd  = (a1[0]-a2[0])*(a1[0]-a2[0]);
  dd += (a1[1]-a2[1])*(a1[1]-a2[1]);
  dd += (a1[2]-a2[2])*(a1[2]-a2[2]);
  return sqrt(dd);
}

double dist2(int x1, int y1, int z1, int x2, int y2, int z2, double sep, int nx, int ny, int nz, int nbins, int inout[]) {
  int l;
  int lx,ly,lz;
  double dd;
  double pos1[3], pos2[3];
  pos1[0] = sep*(double)x1 + sep/2.0;
  pos1[1] = sep*(double)y1 + sep/2.0;
  pos1[2] = sep*(double)z1 + sep/2.0;
  pos2[0] = sep*(double)x2 + sep/2.0;
  pos2[1] = sep*(double)y2 + sep/2.0;
  pos2[2] = sep*(double)z2 + sep/2.0;
  dd = dist(pos1,pos2);
  for(l=0;l<nbins;l++) {
    lx = (int)(((pos2[0]-pos1[0])*(double)l/((double)nbins-1.0)+pos1[0])-sep/2.0)/sep;
    ly = (int)(((pos2[1]-pos1[1])*(double)l/((double)nbins-1.0)+pos1[1])-sep/2.0)/sep;
    lz = (int)(((pos2[2]-pos1[2])*(double)l/((double)nbins-1.0)+pos1[2])-sep/2.0)/sep;
    if(inout[ny*nz*lx+nz*ly+lz] == 1) {
      dd = 10000.0;
      break;
    }
  }
  return dd;
}

int mindistance(double dist[], int sptSet[], int nbins) {
   double min = 10000.0;
   int min_index;
   for(int v=0;v<nbins;v++) {
     if (sptSet[v] == 0 && dist[v] <= min) {
       min = dist[v];
       min_index = v;
     }
   }
   return min_index;
}

int main (int argc, char *argv[]) {
  int i, j, k, l, m, u;
  int ix, iy, iz;
  int jx, jy, jz;
  int nx, ny, nz;
  int cx, cy, cz;
  int natoms, inside, maxn;
  int ben[3];
  int inout[100000];
  int sptSet[100000];
  FILE* fp;
  char elems[5] = {'H','C','N','O','S'};
  char buf[200];
  char subbuf[3];
  char subbuf2[8];
  char atoms[100000];
  double d;
  double sep;
  double vdws[5];
  double coords[100000][3];
  double pos[3];
  double pos2[3];
  double vdw2;
  double distance[100000];

  vdws[0] = 1.2;
  vdws[1] = 1.7;
  vdws[2] = 1.55;
  vdws[3] = 1.52;
  vdws[4] = 1.8;

  if(argc!=9) {
    printf("Error: wring number of input parameters.\n");
    printf("Usage: ./distancefield.o protein.pdb 70 70 70 2 10 10 10 > distances.txt\n");
    exit(0);
  }

  sscanf(argv[2],"%d",&nx);
  sscanf(argv[3],"%d",&ny);
  sscanf(argv[4],"%d",&nz);
  sscanf(argv[5],"%lf",&sep);
  sscanf(argv[6],"%d",&cx);
  sscanf(argv[7],"%d",&cy);
  sscanf(argv[8],"%d",&cz);

  maxn = nx;
  if(ny>maxn) maxn = ny;
  if(nz>maxn) maxn = nz;

  fp = fopen(argv[1], "r");
  i = 0;
  while (fgets(buf, sizeof(buf), fp) != NULL) {
    buf[strlen(buf) - 1] = '\0';
    if(!strcmp(memcpy(subbuf,&buf[0],4),"ATOM")) {
      atoms[i] = buf[13];
      sscanf(memcpy(subbuf2,&buf[30],8), "%lf", &d);
      coords[i][0] = d;
      sscanf(memcpy(subbuf2,&buf[38],8), "%lf", &d);
      coords[i][1] = d;
      sscanf(memcpy(subbuf2,&buf[46],8), "%lf", &d);
      coords[i][2] = d;
      i++;
    }
  }
  fclose(fp);
  natoms = i;

  for(i=0;i<nx*ny*nz;i++) {
    iz = i % nz;
    iy = i % (ny*nz)/nz;
    ix = i % (nx*ny*nz)/ny/nz;
    pos[0] = sep*(double)ix + sep/2.0;
    pos[1] = sep*(double)iy + sep/2.0;
    pos[2] = sep*(double)iz + sep/2.0;
    inside = 0;
    for(l=0;l<natoms;l++) {
      for(m=0;m<5;m++) {
        if(atoms[l]==elems[m]) {
          vdw2 = vdws[m]+1.4;
        }
      }
      if(dist(pos, coords[l]) < vdw2) {
        inside = 1;
        break;
      }
    }
    inout[i] = inside;
  }

/*
  // PRINT grid.pdb
  for(i=0;i<nx;i++) {
    for(j=0;j<ny;j++) {
      for(k=0;k<nz;k++) {
        if(inout[ny*nz*i+nz*j+k] == 1) {
          pos[0] = 2.0*(double)i + 1.0;
          pos[1] = 2.0*(double)j + 1.0;
          pos[2] = 2.0*(double)k + 1.0;
          printf("ATOM      1  X1  XXX     1     %7.3f %7.3f %7.3f\n",pos[0],pos[1],pos[2]);
        }
      }
    }
  }
*/

  for(i=0;i<nx*ny*nz;i++) {
    distance[i] = 10000.0;
    sptSet[i] = 0;
  }
  distance[ny*nz*cx+nz*cy+cz] = 0.0;

  for(i=0;i<nx*ny*nz-1;i++) {
    u = mindistance(distance, sptSet, nx*ny*nz);
    iz = u % nz;
    iy = u % (ny*nz)/nz;
    ix = u % (nx*ny*nz)/ny/nz;
    sptSet[u] = 1;
    for(j=0;j<nx*ny*nz;j++) {
      jz = j % nz;
      jy = j % (ny*nz)/nz;
      jx = j % (nx*ny*nz)/ny/nz;
      if (sptSet[j]==0 && dist2(jx, jy, jz, ix, iy, iz, sep, nx, ny, nz, maxn, inout) < 10000.0 && distance[u] != 10000.0 && distance[u]+dist2(jx, jy, jz, ix, iy, iz, sep, nx, ny, nz, maxn, inout) < distance[j]) {
        distance[j] = distance[u] + dist2(ix, iy, iz, jx, jy, jz, sep, nx, ny, nz, maxn, inout);
      }
    }
  }

  for(i=0;i<nx;i++) {
    for(j=0;j<ny;j++) {
      for(k=0;k<nz;k++) {
        if(inout[ny*nz*i+nz*j+k] == 0) {
          printf("%f %f %f %f\n", sep*(double)i+sep/2.0, sep*(double)j+sep/2.0, sep*(double)k+sep/2.0, distance[ny*nz*i+nz*j+k]);
        }
      }
    }
  }

  return 0;
}


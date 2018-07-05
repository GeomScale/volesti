// Copyright(c) 2016-2018 Kjell Konis <kjell.konis@me.com>.
// Version: 5.5.2.0-17
// Description: The lpSolveAPI package provides an R interface to 'lp_solve',
// a Mixed Integer Linear Programming (MILP) solver with support for pure
//        linear, (mixed) integer/binary, semi-continuous and special ordered sets
//        (SOS) models.
// License: LGPL-2
// Repository: CRAN

#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "lp_lib.h"

#include "ini.h"

FILE *ini_create(char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");

  return(fp);
}

FILE *ini_open(char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");

  return(fp);
}

void ini_writecomment(FILE *fp, char *comment)
{
  fprintf(fp, "; %s\n", comment);
}

void ini_writeheader(FILE *fp, char *header, int addnewline)
{
  if((addnewline) && (ftell(fp) > 0))
    fputs("\n", fp);
  fprintf(fp, "[%s]\n", header);
}

void ini_writedata(FILE *fp, char *name, char *data)
{
  if(name != NULL)
    fprintf(fp, "%s=%s\n", name, data);
  else
    fprintf(fp, "%s\n", data);
}

int ini_readdata(FILE *fp, char *data, int szdata, int withcomment)
{
  int l;
  char *ptr;

  if(fgets(data, szdata, fp) == NULL)
    return(0);

  if(!withcomment) {
    ptr = strchr(data, ';');
    if(ptr != NULL)
      *ptr = 0;
  }

  l = (int) strlen(data);
  while((l > 0) && (isspace(data[l - 1])))
    l--;
  data[l] = 0;
  if((l >= 2) && (data[0] == '[') && (data[l - 1] == ']')) {
    memcpy(data, data + 1, l - 2);
    data[l - 2] = 0;
    return(1);
  }
  return(2);
}

void ini_close(FILE *fp)
{
  fclose(fp);
}

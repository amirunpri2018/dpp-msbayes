/*
 * rand2new.c
 *
 * Copyright (C) 2006  Michael Hickerson
 *
 * This file is a part of msDQH, distributed with msBayes.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */
/* modified by Naoki Takebayashi */

/*  Link in this file for random number generation with rand() */
/*changed rand to random 11/4/04*/
#define RAN random()/((double)RAND_MAX+1)	// Mac OS X or BSD Unix (Linux)

#include <stdio.h>
#include <stdlib.h>
#define randomfunction() RAN


double
ran1 ()
{
  return (randomfunction ());
}


unsigned int
seedit (const char *flag, long seed)
{

}

/* unsigned int */
/* seedit (const char *flag) */
/* { */
/*   FILE *fopen (), *pfseed; */
/*   unsigned int seed2; */
/*   if (flag[0] == 't') */
/*     { */
/*       pfseed = fopen ("seedms", "r"); */
/*       if (pfseed == NULL) */
/* 	{ */
/* 	  seed2 = 59243; */
/* 	} */
/*       else */
/* 	{ */
/* 	  fscanf (pfseed, " %d", &seed2); */
/* 	  fclose (pfseed); */
/* 	} */
/*       srandom (seed2); */
/*       return (seed2); */
/*     } */
/*   else */
/*     { */
/*       pfseed = fopen ("seedms", "w"); */
/*       fprintf (pfseed, "%ld \n", random ()); */
/*       return (0); */
/*     } */
/* } */

/** @file configuration.h
 *  @brief Function prototypes for configuration functions.
 */
/* ----------------------------------------------------------------------------
   libconfig - A library for processing structured configuration files
   Copyright (C) 2005-2010  Mark A Lindner

   This file is part of libconfig.

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   as published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, see
   <http://www.gnu.org/licenses/>.
   ----------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <libconfig.h>

#ifndef _CSTREAM_H
# include "cstream.h"
#endif

#ifndef _PARAMETERS_H_
# include "parameters.h"
#endif

#ifndef _SPECIES_H
# include "species.h"
#endif

#ifndef _REACTION_H
# include "reaction.h"
#endif

#define COMPARE_LIMIT 100
#define NO_PARAMETERS 100

typedef struct array_name {
  const char *name;
} ARRAY_NAME ;

ARRAY_NAME parameter_names[NO_PARAMETERS];

bool
change_cfg_integer(const char* name,const char *type,const char *comment,
                   const char *sstring, int *par,int value,int count,
                   bool *change, config_setting_t *setting_default,int i);
bool
change_cfg_double(const char* name,const char *type,const char *comment,
                  const char *sstring,double *par,double value,int count,
                  bool *change, config_setting_t *setting_default,int i);
//bool
//change_cfg_bool(const char* name,const char *type,const char *comment,
//                const char *sstring, bool *par,bool value,int count,
//                bool *change, config_setting_t *setting_default,int i);
bool
change_cfg_string(const char *name,const char *type,const char *comment,
                  const char *sstring,char** par,char* value,int count,
                  bool *change, config_setting_t *setting_default,int i);
void
change_cfg_parameters (const char *name,const char *type,const char *comment,
                       int ii,int ivalue, double dvalue, bool bbool,
                       const char* astring, config_setting_t *setting,
                       int count,bool *change);
bool
read_parameter(config_setting_t *setting1,config_setting_t *setting2,
               int ii,int count,bool *change);

bool
read_specie(config_setting_t *setting,int ii,species_t *temp_s);

bool
read_seed(config_setting_t *setting,int ii,seed_t *temp_se);

bool
read_reaction(config_setting_t *setting,int ii,reaction_t *temp_r);


/** @file rt.c
 *  @brief Module of the Kinetics. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <libconfig.h>

#include "parameters.h"
#include "proto.h"
#include "species.h"
#include "reaction.h"

#define MAX_SPECIES 15
#define MAX_REACTIONS 30
#define MAX_SEEDS 15

/** @brief Photoionization */
reaction_t react_photo = {TRUE, 0, 0, {0}, {0}, {""}, {""}, NULL, NULL, NULL, 0.0, NULL};

species_t  *spec_index[MAX_SPECIES];
reaction_t *reaction_index[MAX_REACTIONS];
seed_t     *seed_index[MAX_SEEDS];
int        no_reactions;

/** @brief Any initialization required by the kinetic model, has to be done
 *  in this function.
 */
void
kinetic_init (void)
{
  int cnt;
  const char *filename;

  read_input_file(kin_input,filename);
  
  if (has_photoionization == 1) 
  {
    react_add(&react_photo);
    ions = find_species_by_name("dummyplus");
    photo_ions = find_species_by_name("xplus");
   }

  printf("\n");
  for (cnt = no_reactions; cnt > 0; cnt--)
      react_add(reaction_index[cnt-1]);

  electrons = find_species_by_name("electrons");
}

/** @brief Reads the kinetic parameters (species, reactions) from file
 *  @a filename */
void
read_input_file(const char *f_kinetic_name, const char *filename)
{
  config_t         cfg_kinetic;
  config_setting_t *setting;
  species_t        *temp_s;
  reaction_t       *temp_r;
  seed_t           *temp_se;
  int              i,cnt,cnt2;

  config_init(&cfg_kinetic);

  /* Read the file f_kinetic_name. If there is an error, report it and exit. */
  if(! config_read_file(&cfg_kinetic,f_kinetic_name))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg_kinetic),
            config_error_line(&cfg_kinetic), config_error_text(&cfg_kinetic));
    config_destroy(&cfg_kinetic);
    exit(EXIT_FAILURE);
  }

  /* Get the program name. */
  if(config_lookup_string(&cfg_kinetic, "filename", &filename))
    printf("Example : %s\n\n", filename);
  else
    fprintf(stdout, "No ' program name' setting in configuration file.\n");

    /* Output a list of all species parameters in the file. */
  setting = config_lookup(&cfg_kinetic, "species");
  if(setting != NULL)
  {
    no_species = config_setting_length(setting);

    printf("# species   = %i\n",no_species);

    for(i = 0; i < no_species; ++i)
    {
      temp_s = (species_t*) malloc(sizeof(species_t));
      temp_s->charge = 0.0;
      temp_s->mass   = 0.0;
      temp_s->name   = "";
      if (read_specie(setting,i,temp_s))
              spec_index[i]=temp_s;
      printf("Species '%10s' has mass %10.1e and charge %10.1e\n",
		      spec_index[i]->name,
		      spec_index[i]->mass,
		      spec_index[i]->charge);
      continue;
    }
  }

  /* Output a list of all seed parameters in the file. */
  setting = config_lookup(&cfg_kinetic, "seed");
  if(setting != NULL)
  {
    no_seed = config_setting_length(setting);

    printf("# seed      = %i\n",no_seed);
    printf("\n");

    for(i = 0; i < no_seed; ++i)
    {
      temp_se = (seed_t*) malloc(sizeof(seed_t));
      temp_se->species = -1;
      temp_se->value   = 0.0;
      temp_se->type    = -1;
      temp_se->x0      = 0.0;
      temp_se->y0      = 0.0;
      temp_se->z0      = 0.0;
      temp_se->sigma_x = 0.0;
      temp_se->sigma_y = 0.0;
      temp_se->sigma_z = 0.0;
      if (read_seed(setting,i,temp_se))
	      seed_index[i]=temp_se;
      printf("Found a seed of species %s and type %s\n",
		      seed_index[i]->kind_species,
		      seed_index[i]->kind_type);
      printf("It has value %10.1e, z-position %10.1e and y-sigma %10.1e\n",
		      seed_index[i]->value,
		      seed_index[i]->z0,
		      seed_index[i]->sigma_y);
      continue;
    }
  }

  /* Output a list of all reaction parameters in the file. */
  setting = config_lookup(&cfg_kinetic, "reactions");
  if(setting != NULL)
  {
    no_reactions = config_setting_length(setting);

    printf("# reactions = %i\n",no_reactions);

    for(i = 0; i < no_reactions; ++i)
    {
      temp_r = (reaction_t*) malloc(sizeof(reaction_t));
      temp_r->is_photo = 0;
      temp_r->nin = 0;
      temp_r->nout = 0;
      temp_r->f = NULL;
      temp_r->rt = NULL;
      temp_r->tablefile = "";
      temp_r->k = 0.0;
      temp_r->next = NULL;
      if (read_reaction(setting,i,temp_r))
        reaction_index[i]=temp_r;

      continue;
    }
  }

  printf("\n");
  for (cnt = 0; cnt < no_reactions; cnt++)
  {
    printf("Reaction #%i has %i input-species, %i output-species.\n",
		    cnt,reaction_index[cnt]->nin,reaction_index[cnt]->nout);
    printf("The rates are defined in %s\n",
		    reaction_index[cnt]->tablefile);
    printf("Input species are:\t");
    for (cnt2 = 0; cnt2 < reaction_index[cnt]->nin; cnt2++) {
      printf("(%i) %s\t",
		      reaction_index[cnt]->input[cnt2],
		      reaction_index[cnt]->inname[cnt2]);
    }
    printf("\nOutput species are:\t");
    for (cnt2 = 0; cnt2 < reaction_index[cnt]->nout; cnt2++) {
      printf("(%i) %s\t",reaction_index[cnt]->output[cnt2],
		         reaction_index[cnt]->outname[cnt2]);
    }
    printf("\n");
  }
}

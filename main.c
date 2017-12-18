#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "io.h"
#include "particles.h"

#define XASSERT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            fflush(stdout);                                             \
            exit(EXIT_FAILURE);                                         \
        } \
    } while (0)

#define	 CUBE(x) (x*x*x)

// Local Variables //

char *fin_ids, *fin_snapshot, *foutbase;
int32_t num_snapshot_files, num_match_files;

#ifdef MPI
int ThisTask, NTask, nodeNameLen;
char *ThisNode;
#endif

// Proto-types //

int32_t parse_params(int32_t argc, char **argv);
int32_t init();
int32_t read_match_header(char *fname, int32_t *N_files, int32_t *i_file, int64_t *ids_ThisFile, int64_t *ids_Total, int32_t *ids_Total_HighWord);
int32_t get_num_match_files(char *fin_ids, int32_t *num_match_files);

// Functions //

void bye()
{
#ifdef MPI
  MPI_Finalize();
  free(ThisNode);
#endif
}

int cmpfunc (const void * a, const void * b) {
   return ( *(int64_t*)a - *(int64_t*)b );
}

int32_t parse_params(int32_t argc, char **argv)
{

  if (argc != 4)
  {
    fprintf(stderr, "Usage : ./particle_linker <Input Particle ID Base> <Input Snapshot Base> <Output Base>\n"); 
    return EXIT_FAILURE; 
  }

  fin_ids = strdup(argv[1]);
  fin_snapshot = strdup(argv[2]);
  foutbase = strdup(argv[3]); 

  printf("==================================================\n");
  printf("Running with parameters :\nIDs to Match : %s\nInput Snapshot : %s\nOutput Path : %s\nn", fin_ids, fin_snapshot, foutbase);
  printf("==================================================\n\n");

  return EXIT_SUCCESS;
}

int32_t init()
{

  int32_t status;
  char fname[1024];
  hdf5_header file_header;

  snprintf(fname, 1024, "%s.0.hdf5", fin_snapshot);
  if (access(fname, F_OK) == -1)
  { 
    fprintf(stderr, "Could not find file %s\n", fname);
    return EXIT_FAILURE;
  }

  file_header = malloc(sizeof(struct hdf5_header_struct));
  if (file_header == NULL)
  { 
    fprintf(stderr, "Could not allocate memory for the snapshot header\n");
    return EXIT_FAILURE;
  }

  status = get_header_params(fname, file_header);
  if (status == EXIT_FAILURE)
  { 
    return EXIT_FAILURE;
  }

  num_snapshot_files = file_header->NumFilesPerSnapshot;
  free(file_header);

  printf("Snapshot particles are split up over %d files\n", num_snapshot_files);

  status = get_num_match_files(fin_ids, &num_match_files);
  if (status == EXIT_FAILURE)
  {
   return EXIT_FAILURE; 
  }

  printf("The match IDs are split up over %d files\n", num_match_files);

  return EXIT_SUCCESS;

}

int32_t read_match_header(char *fname, int32_t *i_file, int32_t *N_files, int64_t *ids_ThisFile, int64_t *ids_Total, int32_t *ids_Total_HighWord)
{

  FILE *header_file;

  header_file = fopen(fname, "rb");
  if (header_file == NULL)
  { 
    fprintf(stderr, "Could not open %s for reading\n", fname);
    return EXIT_FAILURE;
  } 

  fread(i_file, sizeof(int32_t), 1, header_file); 
  fread(N_files, sizeof(int32_t), 1, header_file); 
  fread(ids_ThisFile, sizeof(int64_t), 1, header_file); 
  fread(ids_Total, sizeof(int64_t), 1, header_file); 
  fread(ids_Total_HighWord, sizeof(int32_t), 1, header_file); 

  fclose(header_file);

  return EXIT_SUCCESS;

}

int32_t get_num_match_files(char *fin_base, int32_t *num_match_files)
{

  char fname[1024];
  int32_t i_file, ids_Total_HighWord;
  int64_t ids_ThisFile, ids_Total; 

  int32_t status;

  snprintf(fname, 1024, "%s.0", fin_ids);
  if (access(fname, F_OK) == -1)
  {
    fprintf(stderr, "Could not find file %s\n", fname);
    return EXIT_FAILURE;
  }
  
  status = read_match_header(fname, &i_file, num_match_files, &ids_ThisFile, &ids_Total, &ids_Total_HighWord); 
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;

}

int32_t read_match_ids(char *fin_base, int32_t i_file, int64_t *match_ids)
{

  char fname[1024];
  FILE *match_file;
  
  int32_t N_files, ids_Total_HighWord;
  int64_t ids_ThisFile, ids_Total; 

  snprintf(fname, 1024, "%s.%d", fin_base, i_file);
  if (access(fname, F_OK) == -1)
  {
    fprintf(stderr, "Could not find file %s\n", fname);
    return EXIT_FAILURE;
  }

  match_file = fopen(fname, "rb");
  if (match_file == NULL)
  { 
    fprintf(stderr, "Could not open %s for reading\n", fname);
    return EXIT_FAILURE;
  } 

  // Read the header //

  fread(&i_file, sizeof(int32_t), 1, match_file); 
  fread(&N_files, sizeof(int32_t), 1, match_file); 
  fread(&ids_ThisFile, sizeof(int64_t), 1, match_file); 
  fread(&ids_Total, sizeof(int64_t), 1, match_file); 
  fread(&ids_Total_HighWord, sizeof(int32_t), 1, match_file); 

  // Now read the IDs //

  fread(match_ids, sizeof(int64_t), ids_ThisFile, match_file); 
  
  fclose(match_file);

  return EXIT_SUCCESS;

}

void match_particles(int64_t *match_ids, part_t snapshot_particles_local, part_t matched_particles)
{
printf("HELLO\n");
}

int main(int argc, char **argv)
{

  char buf[1024];
  int32_t status, snapshot_file_idx, match_file_idx;

  int32_t N_match_ids_highword, header_match_file_idx;
  int64_t N_match_ids_thisfile, N_match_ids_total;
  int64_t *match_ids;
  part_t matched_particles;

  int32_t dummy_short;
 
  int32_t i;
 
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  ThisNode = malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));

  MPI_Get_processor_name(ThisNode, &nodeNameLen);
  if (nodeNameLen >= MPI_MAX_PROCESSOR_NAME)
  {
    printf("Node name string not long enough!...\n");
    exit(EXIT_FAILURE);
  }
#endif

  atexit(bye);
 
  status = parse_params(argc, argv); // Set the input parameters.
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }

  status = init();
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }

  num_snapshot_files = 1;  
#ifdef MPI
  for (snapshot_file_idx = ThisTask; snapshot_file_idx < num_snapshot_files; snapshot_file_idx += NTask)
#else
  for (snapshot_file_idx = 0; snapshot_file_idx < num_snapshot_files; ++snapshot_file_idx)
#endif 
  {  
    part_t snapshot_particles_local;
    snapshot_particles_local = malloc(sizeof(struct particle_struct));

    if (snapshot_particles_local == NULL)
    {
      fprintf(stderr, "Could not allocate memory for local snapshot particles for file number %d\n", snapshot_file_idx);
      exit(EXIT_FAILURE);
    }

    fill_particles(fin_snapshot, snapshot_file_idx, snapshot_particles_local);

    // As a 'worst' case scenario we will assume that ALL particles in the snapshot will be matched. // 
    matched_particles = malloc(sizeof(struct particle_struct));

    if (matched_particles == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the matched particles for file number %d\n", snapshot_file_idx);
      exit(EXIT_FAILURE);
    }

    for (i = 0; i < 6; ++i)
    {
      matched_particles->NumParticles_Total[i] = snapshot_particles_local->NumParticles_Total[i];
    }
    matched_particles->NumParticles_Total_AllType = snapshot_particles_local->NumParticles_Total_AllType;
    allocate_particle_memory(matched_particles, matched_particles->NumParticles_Total_AllType);
 
    // Now we have the particles for the snapshot allocated, let's read in the match IDs and start matching! //
        
    num_match_files = 1; 
    for (match_file_idx = 0; match_file_idx < num_match_files; ++match_file_idx)
    { 
      snprintf(buf, 1024, "%s.%d", fin_ids, match_file_idx);
      if (access(buf, F_OK) == -1)
      {
        fprintf(stderr, "Could not find file %s\n", buf);
        exit(EXIT_FAILURE); 
      }

      read_match_header(buf, &header_match_file_idx, &dummy_short, &N_match_ids_thisfile, &N_match_ids_total, &N_match_ids_highword); 
      XASSERT(header_match_file_idx == match_file_idx, "The header for the match ID file said this is file %d but the loop says it is file %d\n", header_match_file_idx, match_file_idx); 

      match_ids = malloc(sizeof(int64_t) * N_match_ids_thisfile);
      if (match_ids == NULL)
      {
        fprintf(stderr, "Could not allocate memory for the match IDs of file %d\n", match_file_idx);
        exit(EXIT_FAILURE);
      }

      read_match_ids(fin_ids, match_file_idx, match_ids);
      printf("Read in %ld match IDs now sorting them\n", N_match_ids_thisfile);
      qsort(match_ids, N_match_ids_thisfile, sizeof(int64_t), cmpfunc);

      match_particles(match_ids, snapshot_particles_local, matched_particles);
      free(match_ids); 

    }

    free_localparticles(&matched_particles);
    free_localparticles(&snapshot_particles_local);
   
  }

  free(fin_ids);
  free(fin_snapshot);
  free(foutbase); 

  // Read in the snapshot particles.
  // Loop over the matching particle IDs
  // Match
  // End Loop; write out shit.

 
  return EXIT_SUCCESS;
} 

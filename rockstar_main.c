#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#include <hdf5.h>
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
int32_t write_hdf5_header(hid_t group_id, int64_t NumPart, int64_t N_match_ids_total, int32_t snapshot_file_idx);
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
  printf("Running with parameters :\nIDs to Match : %s\nInput Snapshot : %s\nOutput Path : %s\n", fin_ids, fin_snapshot, foutbase);
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

void match_particles(int64_t *match_ids, int64_t N_match_ids, part_t snapshot_particles, part_t matched_particles, int64_t *offset)
{
  // Here match_ids is the sorted list of particle IDs that we want to find in the snapshot particle struct.
  // Since match_ids is sorted we will iterate over every particle in the snapshot and search for a match in match_ids. 

  int32_t is_found;
  int64_t count, search_idx, snapshot_particle_idx, particle_idx;

  printf("We are performing matching over %ld particles\n", snapshot_particles->NumParticles_Total[1]);

#ifdef OPENMP
  #pragma omp parallel for private(is_found, count, search_idx, particle_idx)
#endif
  // Since the snapshot_particles contain ALL particle types we will only loop over the halo particles (particle type 1).
  for (snapshot_particle_idx = snapshot_particles->NumParticles_Total[0]; snapshot_particle_idx < (snapshot_particles->NumParticles_Total[0] + snapshot_particles->NumParticles_Total[1]); ++snapshot_particle_idx) 
  {
#ifdef SHOW_MATCH_PROGRESS
    if ((snapshot_particle_idx - snapshot_particles->NumParticles_Total[0]) %  (int64_t) (1e6) == 0)
    {
      printf("%.2f%% done\n", (float) (snapshot_particle_idx - snapshot_particles->NumParticles_Total[0]) / (float) (snapshot_particles->NumParticles_Total[1]) * 100.0); 
    }
#endif 
    is_found = 0; // Tracks when (if) we find the particle.
    count = 0;
    search_idx = ceil(N_match_ids / 2.0); // This is where we start our search in the match IDs array. 
    particle_idx = snapshot_particles->ID[snapshot_particle_idx];  
 
    // We now search through the match IDs for the snapshot particle ID.
    // Since the match IDs is sorted, we start at the middle of the match list. Depending if the match ID is smaller/larger than the snapshot ID, we multiply the search index by 0.5 or 1.5 
    // This method is faster than a simple N^2 search.

    while (is_found == 0) 
    {
      ++count;
      if (particle_idx == match_ids[search_idx]) 
      {
        is_found = 1;
      }
      else if (N_match_ids / pow(2, count) < 1.0) // The smallest index movement is less than 1.  The snapshot ID isn'in the the match list! 
      {
        break;
      } 
      else if (particle_idx > match_ids[search_idx]) // The snapshot ID is larger than the match ID, move down the list. 
      {
        search_idx = search_idx + ceil(N_match_ids / pow(2, count + 1));
      }
      else
      {
        search_idx = search_idx - ceil(N_match_ids / pow(2, count + 1)); // The snapshot ID is smaller than the match ID, move up the list. 
      }

      if (search_idx >= N_match_ids) // Fix edge case.
      {
        search_idx = N_match_ids -1;
      }
     
    }
    if (is_found == 1) // We found the particle so copy its properties over.
    {
      matched_particles->ID[*offset] = snapshot_particles->ID[snapshot_particle_idx]; 
      matched_particles->mass[*offset] = snapshot_particles->mass[snapshot_particle_idx]; 

      matched_particles->posx[*offset] = snapshot_particles->posx[snapshot_particle_idx]; 
      matched_particles->posy[*offset] = snapshot_particles->posy[snapshot_particle_idx]; 
      matched_particles->posz[*offset] = snapshot_particles->posz[snapshot_particle_idx]; 

      matched_particles->vx[*offset] = snapshot_particles->vx[snapshot_particle_idx]; 
      matched_particles->vy[*offset] = snapshot_particles->vy[snapshot_particle_idx]; 
      matched_particles->vz[*offset] = snapshot_particles->vz[snapshot_particle_idx]; 

      ++(*offset);
    }
 
  }

  printf("Found %ld matching particles\n", *offset);
}

int32_t write_matched_particles(part_t matched_particles, int64_t NumPart, int64_t N_match_ids_total, int32_t snapshot_file_idx)
{

  char fname[1024]; 
  int64_t *buffer_array_long, i;
  double *buffer_array_double;
  double **buffer_array_multi;

  hid_t       file_id, dataset_id, dataspace_id, group_id, datatype;  /* identifiers */
  hsize_t     dims[2];
  herr_t      status;

  snprintf(fname, 1024, "%s.%d.hdf5", foutbase, snapshot_file_idx);

  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id == -1)
  {
    fprintf(stderr, "Cannot open file %s for writing\n", fname);
    return EXIT_FAILURE;  
  }

  dims[0] = NumPart; // This is the number of particles we are writing.

  buffer_array_long = malloc(sizeof(int64_t) * NumPart);
  if (buffer_array_long == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for buffer array long\n");
    return EXIT_FAILURE;
  }
 
  buffer_array_double = malloc(sizeof(double) * NumPart);
  if (buffer_array_double == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for buffer array double\n");
    return EXIT_FAILURE;
  }

  buffer_array_multi = (double **)malloc(sizeof(double *) * NumPart); // Malloc the top level array.
  buffer_array_multi[0] = (double *)malloc(sizeof(double) * NumPart*3); // Then allocate a contiguous block of memory for the particles.

  for (i = 1; i < NumPart; ++i) buffer_array_multi[i] = buffer_array_multi[0]+ i*3;

  group_id = H5Gcreate2(file_id, "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
 
  // Let's start with the particle IDs //

  for (i = 0; i < NumPart; ++i)
  {
    buffer_array_long[i] = matched_particles->ID[i];
//    printf("%ld\n", buffer_array_long[i]);
  } 

  dataspace_id = H5Screate_simple(1, &dims[0], NULL); // Creates the dataspace.

  dataset_id = H5Dcreate2(file_id, "/PartType1/ParticleIDs", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // Create the dataset inside the group. 
  status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer_array_long); // Then write the particle IDs.
  status = H5Dclose(dataset_id); // End access to the dataset and release resources used by it.
  status = H5Sclose(dataspace_id); // Terminate access to the data space.

  printf("Successfully wrote the particle IDs.\n");

  // Position //

  for (i = 0; i < NumPart; ++i)
  {
    buffer_array_multi[i][0] = matched_particles->posx[i];
    buffer_array_multi[i][1] = matched_particles->posy[i];
    buffer_array_multi[i][2] = matched_particles->posz[i];

//    printf("x = %.4f \t y = %.4f \t z = %.4f\n", buffer_array_multi[i][0], buffer_array_multi[i][1], buffer_array_multi[i][2]);

  } 

  dims[1] = 3; 
  dataspace_id = H5Screate_simple(2, dims, NULL); 

  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  status = H5Tset_order(datatype, H5T_ORDER_LE); 

  dataset_id = H5Dcreate2(file_id, "/PartType1/Coordinates", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer_array_multi[0][0]); 
  status = H5Dclose(dataset_id); 
  status = H5Sclose(dataspace_id); 

  printf("Successfully wrote the particle positions.\n");

  // Velocity //
 
  for (i = 0; i < NumPart; ++i)
  {
    buffer_array_multi[i][0] = matched_particles->vx[i];
    buffer_array_multi[i][1] = matched_particles->vy[i];
    buffer_array_multi[i][2] = matched_particles->vz[i];

  //  printf("x = %.4f \t y = %.4f \t z = %.4f\n", buffer_array_multi[i][0], buffer_array_multi[i][1], buffer_array_multi[i][2]);

  } 

  dims[1] = 3; 
  dataspace_id = H5Screate_simple(2, dims, NULL); // Creates the dataspace.

  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  status = H5Tset_order(datatype, H5T_ORDER_LE); 

  dataset_id = H5Dcreate2(file_id, "/PartType1/Velocities", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer_array_multi[0][0]); 
  status = H5Dclose(dataset_id); 
  status = H5Sclose(dataspace_id); 

  status = H5Gclose(group_id);

  // All the data written, create some attributes //

  group_id = H5Gcreate(file_id, "/Header", 0, H5P_DEFAULT, H5P_DEFAULT);

  status = write_hdf5_header(group_id, NumPart, N_match_ids_total, snapshot_file_idx);
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  } 
  
  H5Gclose(group_id);
  // Better free everything!

  status = H5Fclose(file_id);

  free(buffer_array_multi[0]);
  free(buffer_array_multi);
  free(buffer_array_double);
  free(buffer_array_long);

  return EXIT_SUCCESS; 

}

int32_t write_hdf5_header(hid_t group_id, int64_t NumPart, int64_t N_match_ids_total, int32_t snapshot_file_idx) 
{
  // Notice here I'm doing something a little bit tricky.
  // Normally we could only write out the local information as the attributes and then at the end of the program I'd need to flick back through all the files and write out the global, 'NumPart_Total' information after gathering the numbers from all processors.
  // However we KNOW that the number of IDs that we MUST match is given by N_match_ids_total; the IDs in the match list MUST be in the particle snapshot.
  // Hence we know the number of particles that will be in written for each global snapshot.

  hsize_t header_dim[1] = { 6 };
  int32_t NumPart_ThisFile[6], NumPart_Total[6], NumPart_Total_HighWord[6], i, status;
  char fname[1024];  
  hdf5_header file_header;

  hid_t dataspace_id, attribute_id;

  // Set up the header arrays for the number of particles in this file // 

  file_header = malloc(sizeof(struct hdf5_header_struct));
  if (file_header == NULL)
  { 
    fprintf(stderr, "Could not allocate memory for the snapshot header\n");
    return EXIT_FAILURE;
  }

  snprintf(fname, 1024, "%s.%d.hdf5", fin_snapshot, snapshot_file_idx); 
  status = get_header_params(fname, file_header);
  if (status == EXIT_FAILURE)
  { 
    return EXIT_FAILURE;
  }
 
  for (i = 0; i < 6; ++i)
  {
    if (i == 1)
    {
      NumPart_ThisFile[i] = NumPart;
      NumPart_Total[i] = N_match_ids_total; 
      NumPart_Total_HighWord[i] = N_match_ids_total >> 32; 
    }
    else
    {
      NumPart_ThisFile[i] = 0;
      NumPart_Total[i] = 0;
      NumPart_Total_HighWord[i] = 0;
    }
  }

  // Time to create and write arrays //

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "NumFilesPerSnapshot", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_INT, &num_snapshot_files);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "BoxSize", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &file_header->BoxSize);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "Time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &file_header->Time);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "Redshift", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &file_header->Redshift);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "Omega0", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &file_header->Omega0);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "OmegaLambda", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &file_header->OmegaLambda);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SCALAR);
  attribute_id = H5Acreate(group_id, "HubbleParam", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &file_header->HubbleParam);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace_id, 1, header_dim, NULL);
  attribute_id = H5Acreate(group_id, "NumPart_ThisFile", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_INT, &NumPart_ThisFile); 
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace_id, 1, header_dim, NULL);
  attribute_id = H5Acreate(group_id, "MassTable", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &file_header->MassTable); 
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace_id, 1, header_dim, NULL);
  attribute_id = H5Acreate(group_id, "NumPart_Total", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_INT, &NumPart_Total);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);

  dataspace_id = H5Screate(H5S_SIMPLE);
  H5Sset_extent_simple(dataspace_id, 1, header_dim, NULL);
  attribute_id = H5Acreate(group_id, "NumPart_Total_HighWord", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_INT, &NumPart_Total_HighWord);
  H5Aclose(attribute_id);
  H5Sclose(dataspace_id);


  free(file_header);

  return EXIT_SUCCESS;

}

int32_t do_final_check(int64_t N_matched_total, int64_t N_match_ids_total)
{
  // Here N_matched_total is the number of particles that have been successfully matched by this rank.
  // N_match_ids_total is the total number of particles from the particle match ID files that we were required to match.

  int64_t N_matched_global;

#ifdef MPI  
  MPI_Reduce(&N_matched_total, &N_matched_global, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  N_matched_global = N_matched_total;
#endif

#ifdef MPI
  if (ThisTask == 0)
#endif
  {
    printf("==================================================\n");
    printf("Summary Stats\nNumber particles required to match : %ld\nNumber particles we matched : %ld\n", N_match_ids_total, N_matched_global);
    printf("==================================================\n\n");
    if (N_matched_global != N_match_ids_total)
    {
      fprintf(stderr, "We failed to match all the particles we were required to.\nThese particles MUST be in the snapshot somewhere...\n");
      return EXIT_FAILURE;
    }

  }

  return EXIT_SUCCESS;



}

int main(int argc, char **argv)
{

  char buf[1024];
  int32_t status, snapshot_file_idx, match_file_idx;
  
  int64_t N_matched_thisfile, N_matched_total = 0; // These are the number of IDs that have been matched for each snapshot subfile (_thisfile) and for the entire snapshot (_total).
  
  int32_t N_match_ids_highword, header_match_file_idx;
  int64_t N_match_ids_thisfile, N_match_ids_total; // These are the number of IDs that need to be matched for each snapshot subfile (_thisfile) and for the entire snapshot (_total).
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
    N_matched_thisfile = 0; // Offset for slicing in matched properties from each match ID subfile. 
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

      match_particles(match_ids, N_match_ids_thisfile, snapshot_particles_local, matched_particles, &N_matched_thisfile);
      free(match_ids); 

    }

    status = write_matched_particles(matched_particles, N_matched_thisfile, N_match_ids_total, snapshot_file_idx);
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }

    free_localparticles(&matched_particles);
    free_localparticles(&snapshot_particles_local);
  
    N_matched_total += N_matched_thisfile; 
 
  }

  status = do_final_check(N_matched_total, N_match_ids_total);
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }

  XASSERT(N_match_ids_total == N_matched_total, "From the particle ID match file we required to match %ld particles.  However, after searching all %d snapshot subfiles, we only matched %ld IDs.\nThese matched particles MUST be in the snapshot somewhere...\n", N_match_ids_total, num_snapshot_files, N_matched_total);

  free(fin_ids);
  free(fin_snapshot);
  free(foutbase); 

  // Read in the snapshot particles.
  // Loop over the matching particle IDs
  // Match
  // End Loop; write out shit.

 
  return EXIT_SUCCESS;
} 

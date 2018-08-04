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
int32_t num_snapshot_files, num_fof_files;

int32_t ThisTask, NTask;
#ifdef MPI
int32_t nodeNameLen;
char *ThisNode;
#endif

// Proto-types //

int32_t parse_params(int32_t argc, char **argv);
int32_t init();
int32_t read_match_header(char *fname, int32_t *N_files, int32_t *i_file, int64_t *ids_ThisFile, int64_t *ids_Total, int32_t *ids_Total_HighWord);
int32_t write_hdf5_header(hid_t group_id, int64_t N_matched_ThisSnap, int64_t N_matched_AllSnap, int32_t snapshot_file_idx, int64_t N_fof_ids);
int32_t read_fof_ids_ThisTask(char *fin_base, int32_t ThisTask, int32_t NTask, int64_t **fof_ids_ThisTask, int64_t *N_fof_ids_ThisTask, int64_t *N_fof_ids_AllTask);
int32_t determine_N_fof_ids_ThisTask(char *fin_ids, int32_t ThisTask, int32_t NTask, int64_t *N_fof_ids_thistask);
int32_t read_fof_ids(char *fin_base, int32_t i_file, int64_t **fof_ids, int64_t *Nids_ThisFile, int64_t *Nids_AllFile);
int32_t write_matched_particles_ThisSnap(part_t matched_particles, int64_t N_matched_ThisSnap, int64_t N_matched_AllSnap, int32_t snapshot_file_idx, int32_t ThisTask, int64_t N_fof_ids);
void match_particles(int64_t *match_ids, int64_t N_match_ids, part_t snapshot_particles, part_t matched_particles, int64_t *offset, int64_t *N_matched_ThisSnap);

// Functions //

void bye()
{
#ifdef MPI
  MPI_Finalize();
  free(ThisNode);
#endif
}
static __inline int simple_cmp(const void *a, const void *b) {
  const int64_t da = *((const int64_t *) a);
  const int64_t db = *((const int64_t *) b);
  return (da < db) ? -1 : (da == db) ? 0 : 1;
}


int32_t cmpfunc (const void * a, const void * b) {
  if ( *(int64_t*)a - *(int64_t*)b < 0);
    return -1;
  if ( *(int64_t*)a - *(int64_t*)b > 0);
    return -1;
  return 0;
}

int32_t parse_params(int32_t argc, char **argv)
{

  char check_fof[1024], check_snapshot[1024], check_logfile[1024];
  int32_t snapnum, force_calculation;

  if (argc != 6)
  {
    fprintf(stderr, "This is the version of the particle linker that will link the FoF particles for the Kali simulation.\n");
    fprintf(stderr, "Usage : ./kali_linker <FoF IDs Base> <Kali Snapshot Base> <Output Base> <Snapshot Number> <Force Calculation>\n"); 
    fprintf(stderr, "Since we check to see if a calculation has been done before (using the size of the log file as an indicator), <Force Calculation> = 1 will ignore this check and will do the linking regardless.");
    return EXIT_FAILURE; 
  }

  fin_ids = strdup(argv[1]);

  snprintf(check_fof, 1024, "%s.1008", fin_ids);

  if (access(check_fof, F_OK) == -1)
  { 
    printf("The number of FoF ID files is 1008\n");
    num_fof_files = 1008; 
  }
  else
  {
    printf("The number of FoF ID files is 1540\n");
    num_fof_files = 1540; 
  } 
    
  fin_snapshot = strdup(argv[2]);

  snprintf(check_snapshot, 1024, "%s.0.hdf5", fin_snapshot);
  printf("%s\n", check_snapshot); 
  if (access(check_snapshot, F_OK) == -1)
  { 
    printf("There are no snapshot files for this snapshot.\n");
    return EXIT_FAILURE; 
  }


  foutbase = strdup(argv[3]); 
  snapnum = atoi(argv[4]);

  if (snapnum < 0 || snapnum > 98) 
  {
    fprintf(stderr, "The Kali snapshots range from 0 to 98.  You entered %d as the snapnum.\n", snapnum);
    return EXIT_FAILURE;
  }
  
  force_calculation = atoi(argv[5]);
  
  if ((force_calculation < 0) || (force_calculation > 1))
  {
    fprintf(stderr, "Force calculation can only be 0 or 1.  You entered %d\n", force_calculation);
    return EXIT_FAILURE;
  }

  snprintf(check_logfile, 1024, "%s/logfiles/kali_%03d.log", "/lustre/projects/p134_swin/jseiler/kali/pseudo_snapshots/", snapnum); 

  if (access(check_logfile, F_OK) == -1)
  {
  }
  else
  {
    FILE *logfile;
    int64_t file_size;
    logfile = fopen(check_logfile, "r");
    
    fseek(logfile, 0L, SEEK_END); // Move to the end of the file
    file_size = ftell(logfile); // Then count how many bytes we have in it.
    fclose(logfile);
  
    if (file_size > 320.0 * 1024.0 && file_size < 360.0 * 1024.0)
    {
      fprintf(stderr, "The logfile has a size %.4f KB meaning that it has already been matched.\n", file_size / 1024.0);
      if (force_calculation == 1)
      {
        fprintf(stderr, "However we are forcing the calculation to be done anyway.\n");
      }
      else
      {
        return EXIT_FAILURE;
      }
    }
  } 

  printf("==================================================\n");
  printf("Running with parameters :\nFoF IDs : %s\nNumber of FoF IDs Subfiles : %d\nInput Snapshot : %s\nOutput Path : %s\nSnapshot Number : %d\n", fin_ids, num_fof_files, fin_snapshot, foutbase, snapnum);
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

void match_particles(int64_t *match_ids, int64_t N_fof_ids, part_t snapshot_particles, part_t matched_particles, int64_t *N_matched_AllSnap, int64_t *N_matched_ThisSnap)
{
  // Here match_ids is the sorted list of particle IDs that we want to find in the snapshot particle struct.
  // Since match_ids is sorted we will iterate over every particle in the snapshot and search for a match in match_ids. 

  int32_t is_found;
  int64_t count, search_idx, snapshot_particle_idx, particle_idx, N_matched = 0;
  double progress;

  printf("We are performing matching over %ld particles\n", snapshot_particles->NumParticles_Total[1]);

  int64_t *matched_particle_idx, idx;
  
  matched_particle_idx = malloc(sizeof(int64_t) * N_fof_ids * 1.1); // Assume that at most only a tenth of the FoF particles will be in this snapshot file.
 
  // Since the snapshot_particles contain ALL particle types we will only loop over the halo particles (particle type 1).
  for (snapshot_particle_idx = snapshot_particles->NumParticles_Total[0]; snapshot_particle_idx < (snapshot_particles->NumParticles_Total[0] + snapshot_particles->NumParticles_Total[1]); ++snapshot_particle_idx) 
  {
  
    progress = (snapshot_particle_idx - snapshot_particles->NumParticles_Total[0]) / (snapshot_particles->NumParticles_Total[1]) * 100.0;  
#ifdef SHOW_MATCH_PROGRESS
    if ((snapshot_particle_idx - snapshot_particles->NumParticles_Total[0]) %  (int64_t) (1e6) == 0)
    {
      printf("%.2f%% done\n", progress); 
    }
#endif 
    is_found = 0; // Tracks when (if) we find the particle.
    count = 0;
    search_idx = ceil(N_fof_ids / 2.0); // This is where we start our search in the match IDs array. 
    particle_idx = snapshot_particles->ID[snapshot_particle_idx];  
 
    // We now search through the FoF IDs for the snapshot particle ID.
    // Since the FoF IDs are sorted, we start at the middle of the FoF list. Depending if the FoF ID is smaller/larger than the snapshot ID, we multiply the search index by 0.5 or 1.5 
    // This method is faster than a simple N^2 search.

    /*
    if (particle_idx == 7805847370)
    {
      printf("The particle with ID 7805847370 is now.\n");
      extra_info = 1;
      extra_info = 0;
      
      for (manual_idx = 0; manual_idx < N_fof_ids; ++manual_idx)
      {
        if (match_ids[manual_idx] == 7805847370)
        {
          printf("Manual_idx = %ld\n", manual_idx);

        }
        if (match_ids[manual_idx] > 7805847370)
        {
          printf("Have gone past 7805847370 at manual_idx = %ld. match_ids[manual_idx] = %ld\n", manual_idx, match_ids[manual_idx]);
        }
      }
     
    }
    */
    while (is_found == 0) 
    {

      ++count;

      /*
      if (extra_info == 1)
      {
        printf("%ld\n", match_ids[59189]);
        printf("search_idx = %ld, match_ids[search_idx] = %ld, N_fof_ids = %ld, count = %ld, N_fof_ids / pow(2,count) = %.4f\n", search_idx, match_ids[search_idx], N_fof_ids, count, N_fof_ids / pow(2,count));
      }
      */

      if (particle_idx == match_ids[search_idx]) 
      {
        is_found = 1;
      }
      else if (N_fof_ids / pow(2, count) < 1.0) // The smallest index movement is less than 1.  The snapshot ID isn'in the the match list! 
      {
        break;
      } 
      else if (particle_idx > match_ids[search_idx]) // The snapshot ID is larger than the match ID, move down the list. 
      {
        search_idx = search_idx + ceil(N_fof_ids / pow(2, count + 1));
      }
      else
      {
        search_idx = search_idx - ceil(N_fof_ids / pow(2, count + 1)); // The snapshot ID is smaller than the match ID, move up the list. 
      }

      if (search_idx >= N_fof_ids) // Fix edge case.
      {        
        search_idx = N_fof_ids -1;
      }      
    
      if (search_idx < 0)
      {
        search_idx = 0;
      }
 
    }
        
    if (is_found == 1) // We found the particle so copy its properties over.
    {

      matched_particle_idx[N_matched] = snapshot_particle_idx; 
      ++N_matched;

      XASSERT(N_matched < N_fof_ids * 1.1, "We are %.4f the way through matching the FOF IDs and we have hit the limit of our allocated space for matched_particle_idx.  The limit is %.4f and we are at %ld\n", progress, (int64_t) N_fof_ids * 1.1, N_matched); 
    }

  }

  printf("Found %ld matching particles\n", N_matched);
  printf("Will now allocate memory for the matched particles struct and then grab all the snapshot particle properties.\n");

  allocate_particle_memory(matched_particles, N_matched); 

  for (idx = 0; idx < N_matched; ++idx)
  {
//    printf("idx = %ld \t N_matched = %ld \t snapshot_particle_idx[idx] = %ld\n", idx, N_matched, matched_particle_idx[idx]);
    snapshot_particle_idx = matched_particle_idx[idx]; 

    matched_particles->ID[idx] = snapshot_particles->ID[snapshot_particle_idx]; 
    matched_particles->mass[idx] = snapshot_particles->mass[snapshot_particle_idx]; 

    matched_particles->posx[idx] = snapshot_particles->posx[snapshot_particle_idx]; 
    matched_particles->posy[idx] = snapshot_particles->posy[snapshot_particle_idx]; 
    matched_particles->posz[idx] = snapshot_particles->posz[snapshot_particle_idx]; 

    matched_particles->vx[idx] = snapshot_particles->vx[snapshot_particle_idx]; 
    matched_particles->vy[idx] = snapshot_particles->vy[snapshot_particle_idx]; 
    matched_particles->vz[idx] = snapshot_particles->vz[snapshot_particle_idx]; 
 
  }

  free(matched_particle_idx);

  *N_matched_ThisSnap = N_matched;
  *N_matched_AllSnap += N_matched;

  printf("Everything has been matched!\n");

}

int32_t write_matched_particles_ThisSnap(part_t matched_particles, int64_t N_matched_ThisSnap, int64_t N_matched_AllSnap, int32_t snapshot_file_idx, int32_t ThisTask, int64_t N_fof_ids)
{

  char fname[1024]; 
  int64_t *buffer_array_long, i;
  double *buffer_array_double;
  double **buffer_array_multi;

  hid_t       file_id, dataset_id, dataspace_id, group_id, datatype;  /* identifiers */
  hsize_t     dims[2];
  herr_t      status;

  printf("I am Task %d and I am about to write out the %ld particles I matched this Snapshot. I have currently written out %ld particles.\n", ThisTask, N_matched_ThisSnap, N_matched_AllSnap);

  snprintf(fname, 1024, "%s.%d.hdf5", foutbase, snapshot_file_idx);

  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id == -1)
  {
    fprintf(stderr, "Cannot open file %s for writing\n", fname);
    return EXIT_FAILURE;  
  }

  group_id = H5Gcreate2(file_id, "/PartType1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  dims[0] = N_matched_ThisSnap; // This is the number of particles we are writing.

  if (N_matched_ThisSnap != 0) // Have no IDs matched for this snapshot file so let's write the header and move on.
  {
  
    buffer_array_long = malloc(sizeof(int64_t) * N_matched_ThisSnap);
    if (buffer_array_long == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for buffer array long\n");
      return EXIT_FAILURE;
    }
   
    buffer_array_double = malloc(sizeof(double) * N_matched_ThisSnap);
    if (buffer_array_double == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for buffer array double\n");
      return EXIT_FAILURE;
    }

    buffer_array_multi = (double **)malloc(sizeof(double *) * N_matched_ThisSnap); // Malloc the top level array.
    buffer_array_multi[0] = (double *)malloc(sizeof(double) * N_matched_ThisSnap*3); // Then allocate a contiguous block of memory for the particles.

    for (i = 1; i < N_matched_ThisSnap; ++i) buffer_array_multi[i] = buffer_array_multi[0]+ i*3;
   
    // Let's start with the particle IDs //
    for (i = 0; i < N_matched_ThisSnap; ++i)
    {
      buffer_array_long[i] = matched_particles->ID[i];
    } 

    dataspace_id = H5Screate_simple(1, &dims[0], NULL); // Creates the dataspace.

    dataset_id = H5Dcreate2(file_id, "/PartType1/ParticleIDs", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // Create the dataset inside the group. 
    status = H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer_array_long); // Then write the particle IDs.
    status = H5Dclose(dataset_id); // End access to the dataset and release resources used by it.
    status = H5Sclose(dataspace_id); // Terminate access to the data space.

    // Position //

    for (i = 0; i < N_matched_ThisSnap; ++i)
    {
      buffer_array_multi[i][0] = matched_particles->posx[i];
      buffer_array_multi[i][1] = matched_particles->posy[i];
      buffer_array_multi[i][2] = matched_particles->posz[i];
    } 

    dims[1] = 3; 
    dataspace_id = H5Screate_simple(2, dims, NULL); 

    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE); 

    dataset_id = H5Dcreate2(file_id, "/PartType1/Coordinates", datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer_array_multi[0][0]); 
    status = H5Dclose(dataset_id); 
    status = H5Sclose(dataspace_id); 

    // Velocity //
   
    for (i = 0; i < N_matched_ThisSnap; ++i)
    {
      buffer_array_multi[i][0] = matched_particles->vx[i];
      buffer_array_multi[i][1] = matched_particles->vy[i];
      buffer_array_multi[i][2] = matched_particles->vz[i];

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

  }

  // All the data written, create some attributes //

  group_id = H5Gcreate(file_id, "/Header", 0, H5P_DEFAULT, H5P_DEFAULT);

  status = write_hdf5_header(group_id, N_matched_ThisSnap, N_matched_AllSnap, snapshot_file_idx, N_fof_ids);

  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  } 
  
  H5Gclose(group_id);

  status = H5Fclose(file_id);

  if (N_matched_ThisSnap != 0) // Have no IDs matched for this snapshot file so let's write the header and move on.
  {
    free(buffer_array_multi[0]);
    free(buffer_array_multi);
    free(buffer_array_double);
    free(buffer_array_long);
  }

  return EXIT_SUCCESS; 

}

int32_t write_hdf5_header(hid_t group_id, int64_t N_matched_ThisSnap, int64_t N_matched_AllSnap, int32_t snapshot_file_idx, int64_t N_fof_ids) 
{

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
      NumPart_ThisFile[i] = N_matched_ThisSnap;
      NumPart_Total[i] = N_fof_ids; 
      NumPart_Total_HighWord[i] = N_fof_ids >> 32; 
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

int32_t read_fof_ids_ThisTask(char *fin_base, int32_t ThisTask, int32_t NTask, int64_t **fof_ids_ThisTask, int64_t *N_fof_ids_ThisTask, int64_t *N_fof_ids_AllTask)
{

  int32_t fof_file_idx, status;
  int64_t offset = 0, Nids_ThisFile, fof_idx;
  int64_t *tmp_fof_ids_thistask, *tmp_fof_ids_thisfile;

  status = determine_N_fof_ids_ThisTask(fin_base, ThisTask, NTask, N_fof_ids_ThisTask);
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  tmp_fof_ids_thistask = malloc(sizeof(int64_t) * (*N_fof_ids_ThisTask)); 
  if(tmp_fof_ids_thistask == NULL)
  {
    fprintf(stderr, "Could not allocate enough memory to hold the FoF IDs for task %d\n", ThisTask);
    return EXIT_FAILURE;
  }

  // We have now determined how many FoF IDs this task will read and allocated enough memory to hold them. //
  // Now lets loop over each file that this task will read in, read the FoF IDs and then slice it into the global array. // 

  for (fof_file_idx = ThisTask; fof_file_idx < num_fof_files; fof_file_idx += NTask)
  { 

#ifdef DEBUG_FOF_IDS
    if(fof_file_idx % 10 == 0)
    {
      printf("On FoF File %d\n", fof_file_idx);
    }
#endif 
 
    read_fof_ids(fin_base, fof_file_idx, &tmp_fof_ids_thisfile, &Nids_ThisFile, N_fof_ids_AllTask); 

    for (fof_idx = 0; fof_idx < Nids_ThisFile; ++fof_idx)
    {
    
      tmp_fof_ids_thistask[fof_idx + offset] = tmp_fof_ids_thisfile[fof_idx];

    }
    offset += Nids_ThisFile;
    free(tmp_fof_ids_thisfile);
  } 

  *fof_ids_ThisTask = tmp_fof_ids_thistask;

  return EXIT_SUCCESS;

}

int32_t read_fof_ids(char *fin_base, int32_t i_file, int64_t **fof_ids, int64_t *Nids_ThisFile, int64_t *Nids_AllFile)
{

  char fname[1024];
  FILE *fof_file;
  
  int32_t Ngroups, TotNgroups, Nids, Ntask, Offset; 
  int64_t TotNids, *ids; 

  snprintf(fname, 1024, "%s.%d", fin_base, i_file);
  if (access(fname, F_OK) == -1)
  {
    fprintf(stderr, "Could not find file %s\n", fname);
    return EXIT_FAILURE;
  }

  fof_file = fopen(fname, "rb");
  if (fof_file == NULL)
  { 
    fprintf(stderr, "Could not open %s for reading\n", fname);
    return EXIT_FAILURE;
  } 

  // Read the header //

  fread(&Ngroups, sizeof(int32_t), 1, fof_file); 
  fread(&TotNgroups, sizeof(int32_t), 1, fof_file); 
  fread(&Nids, sizeof(int32_t), 1, fof_file); 
  fread(&TotNids, sizeof(int64_t), 1, fof_file); 
  fread(&Ntask, sizeof(int32_t), 1, fof_file); 
  fread(&Offset, sizeof(int32_t), 1, fof_file); 

  ids = malloc(sizeof(int64_t) * Nids);
  if (ids == NULL)
  {
    fprintf(stderr, "Could not allocate memory for reading of FoF IDs from file %s\n", fname);
    return EXIT_FAILURE;
  } 

  fread(ids, sizeof(int64_t), Nids, fof_file); 

  int64_t i;
  for(i = 0; i < Nids; ++i)
  {
    if (ids[i] == 247991951736)
    {
      printf("%ld\n", ids[i]);
    }
  }
  fclose(fof_file);
 
  *fof_ids = ids; 
  *Nids_ThisFile = (int64_t) Nids;
  *Nids_AllFile = TotNids;

  return EXIT_SUCCESS;


}

int32_t read_fof_ids_header(char *fname, int32_t *Ngroups, int32_t *TotNgroups, int32_t *Nids, int64_t *TotNids, int32_t *Ntask, int32_t *Offset)
{

  FILE *fof_file;

  fof_file = fopen(fname, "rb");
  if (fof_file == NULL)
  {     
    fprintf(stderr, "Could not open %s for reading\n", fname);
    return EXIT_FAILURE;
  } 

  // Read the header //

  fread(Ngroups, sizeof(int32_t), 1, fof_file); 
  fread(TotNgroups, sizeof(int32_t), 1, fof_file); 
  fread(Nids, sizeof(int32_t), 1, fof_file); 
  fread(TotNids, sizeof(int64_t), 1, fof_file); 
  fread(Ntask, sizeof(int32_t), 1, fof_file); 
  fread(Offset, sizeof(int32_t), 1, fof_file); 

  fclose(fof_file);
  return EXIT_SUCCESS;


}

int32_t determine_N_fof_ids_ThisTask(char *fin_ids, int32_t ThisTask, int32_t NTask, int64_t *N_fof_ids_thistask)
{

  char fname[1024];

  int32_t Ngroups, TotNgroups, Nids, Ntask, Offset, fof_file_idx, status; 
  int64_t TotNids;

  *N_fof_ids_thistask = 0;

  printf("I am Task %d and I am determining how many IDs I will be reading in\n", ThisTask);

  for (fof_file_idx = ThisTask; fof_file_idx < num_fof_files; fof_file_idx += NTask)
  { 

    snprintf(fname, 1024, "%s.%d", fin_ids, fof_file_idx);

    status = read_fof_ids_header(fname, &Ngroups, &TotNgroups, &Nids, &TotNids, &Ntask, &Offset);
    if(status == EXIT_FAILURE)
    {
      return EXIT_FAILURE;
    }
     
    *N_fof_ids_thistask += Nids;
  }

  if (NTask == 1 && *N_fof_ids_thistask != TotNids)
  {
    printf("We are only running with 1 Task but we have not read in the all the FoF IDs\nWe read in %ld compared to the total number %ld\n", *N_fof_ids_thistask, TotNids);
    return EXIT_FAILURE;
  }
    
  printf("I will be reading in %ld IDs.  This is %.4f of the total number of FoF IDs.\n", *N_fof_ids_thistask, (double) *N_fof_ids_thistask / TotNids); 
 
  return EXIT_SUCCESS;
}



int32_t main(int argc, char **argv)
{

  int32_t status, snapshot_file_idx, i; 
  int64_t N_matched_ThisSnap = 0, N_matched_AllSnap = 0, N_fof_ids_ThisTask, N_fof_ids_AllTask, ii;
  int64_t *fof_ids; 
  part_t matched_particles;
 
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
#else
  ThisTask = 0;
  NTask = 1;
#endif

  atexit(bye);

  status = parse_params(argc, argv); // Set the input parameters.
  if (status == EXIT_FAILURE)
  {
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(EXIT_FAILURE);
#endif
  }

#ifdef ADJUST_PROCESSORS 
  int32_t max_processors; 

  max_processors = 58; 

  if (ThisTask >= max_processors) 
  {
    fprintf(stderr, "I am Task %d and I am doing nothing.\n", ThisTask);
     
  }
  else
  {
 
    if (NTask > max_processors)
    {
      NTask = max_processors;
    }  
#endif
 
  status = init();
  if (status == EXIT_FAILURE)
  {
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(EXIT_FAILURE);
#endif
  }
 
  // The logic flow of the Kali linker is different to others. //
  // Because there are so many particles in FoFs, a single task cannot hold all of the FoF IDs (and still be below 16gb per task, the limit of largemem). //
  // Instead we will split the FoF IDs over the tasks and then each task will search through each snapshot for their corresponding particles. //
 
  printf("Starting to read in the FoF IDs.\n"); 
  status = read_fof_ids_ThisTask(fin_ids, 0, 1, &fof_ids, &N_fof_ids_ThisTask, &N_fof_ids_AllTask);

  if (status == EXIT_FAILURE)
  {
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    exit(EXIT_FAILURE);
#endif
  }
 
  printf("Read in %ld FoF IDs now sorting them\n", N_fof_ids_ThisTask);   
  qsort(fof_ids, N_fof_ids_ThisTask, sizeof(int64_t), simple_cmp);
  
  for (ii = 0; ii < N_fof_ids_ThisTask-1; ++ii)
  {
    XASSERT(fof_ids[ii] < fof_ids[ii+1], "i = %ld \t fof_ids[i] = %ld \t fof_ids[i+1] = %ld\n", ii, fof_ids[ii], fof_ids[ii+1]);
  } 
  
  // Now lets read in each snapshot and do the matching. //  
  
  for (snapshot_file_idx = ThisTask; snapshot_file_idx < num_snapshot_files; snapshot_file_idx += NTask)
  { 
    part_t snapshot_particles_local;
    snapshot_particles_local = malloc(sizeof(struct particle_struct));
    
    if (snapshot_particles_local == NULL)
    {
      fprintf(stderr, "Could not allocate memory for local snapshot particles for file number %d\n", snapshot_file_idx);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      exit(EXIT_FAILURE);
#endif
    }

    printf("I am Task %d and I am reading from snapshot file %d\n", ThisTask, snapshot_file_idx);
    fill_particles(fin_snapshot, snapshot_file_idx, snapshot_particles_local);

    // Now allocate struct for the matched particles.
    
    matched_particles = malloc(sizeof(struct particle_struct));

    if (matched_particles == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the matched particles\n"); 
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      exit(EXIT_FAILURE);
#endif
    }

    for (i = 0; i < 6; ++i)
    {
      matched_particles->NumParticles_Total[i] = 0; 
    }
    matched_particles->NumParticles_Total[1] = N_fof_ids_ThisTask; // We only have this type of particle for Kali.
    matched_particles->NumParticles_Total_AllType = N_fof_ids_ThisTask; 
  
    match_particles(fof_ids, N_fof_ids_ThisTask, snapshot_particles_local, matched_particles, &N_matched_AllSnap, &N_matched_ThisSnap);
            
    printf("Matched the FoF IDs. Time to write.\n");
    status = write_matched_particles_ThisSnap(matched_particles, N_matched_ThisSnap, N_matched_AllSnap, snapshot_file_idx, ThisTask, N_fof_ids_ThisTask);

    if (status == EXIT_FAILURE)
    {
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      exit(EXIT_FAILURE);
#endif
    }
   
    free_localparticles(&snapshot_particles_local);      
    free_localparticles(&matched_particles);
  }
 
  free(fof_ids); 
  free(fin_ids);
  free(fin_snapshot);
  free(foutbase); 


#ifdef ADJUST_PROCESSORS 
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

/*
  status = do_final_check(N_matched_AllSnap, N_fof_ids_ThisTask);
  if (status == EXIT_FAILURE)
  {
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      exit(EXIT_FAILURE);
#endif
  }
*/
  printf("I am Task %d and I'm all done and leaving!\n", ThisTask); 
  return EXIT_SUCCESS;
} 

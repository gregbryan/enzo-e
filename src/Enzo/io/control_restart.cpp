// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_control_restart.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-10
/// @brief    Enzo-E portion of restart
/// @ingroup  Enzo
///
/// This file controls restarting from files generated by a previous call
/// to EnzoMethodCheck

//--------------------------------------------------

#include "Cello/charm_simulation.hpp"
#include "Cello/charm_mesh.hpp"
#include "Cello/main.hpp"

#include "Enzo/enzo.hpp"
#include "Enzo/io/io.hpp"

//  #define DEBUG_RESTART
//  #define TRACE_BLOCK
// #define PRINT_FIELD_RESTART
//  #define TRACE_SYNC

//--------------------------------------------------
#ifdef TRACE_SYNC
#   undef TRACE_SYNC
#   define TRACE_SYNC(SYNC,MSG)                                         \
  CkPrintf ("TRACE_SYNC %p %s %d/%d\n", \
            (void *)(&SYNC), std::string(MSG).c_str(),SYNC.value(),SYNC.stop()); \
  fflush(stdout);
#else
#   define TRACE_SYNC(SYNC,MSG) /* ... */
#endif
//--------------------------------------------------
#ifdef PRINT_FIELD_RESTART
#   undef PRINT_FIELD_RESTART

#   define PRINT_FIELD_RESTART(MSG,FIELD,DATA)                                  \
  {                                                                     \
    Field field = DATA->field();                                        \
    int mx,my,mz;                                                       \
    int gx,gy,gz;                                                       \
    int index_field = field.field_id(FIELD);                            \
    field.dimensions(index_field,&mx,&my,&mz);                          \
    field.ghost_depth(index_field,&gx,&gy,&gz);                         \
    enzo_float * value = (enzo_float *)field.values(FIELD);             \
    enzo_float min=1e30;                                                \
    enzo_float max=-1e30;                                               \
    enzo_float sum=0;                                                   \
    for (int iz=gz; iz<mz-gz; iz++) {                                   \
      for (int iy=gy; iy<my-gy; iy++) {                                 \
        for (int ix=gx; ix<mx-gx; ix++) {                               \
          const int i=ix+mx*(iy+my*iz);                                 \
          min=std::min(min,value[i]);                                   \
          max=std::max(max,value[i]);                                   \
          sum+=value[i];                                                \
        }                                                               \
      }                                                                 \
    }                                                                   \
    CkPrintf ("PRINT_FIELD_RESTART %s %s  %g %g %g\n",MSG,FIELD,min,max,sum/((mx-2*gx)*(my-2*gy)*(mz-2*gz))); \
    fflush(stdout);                                                     \
  }
#else
#   define PRINT_FIELD_RESTART(MSG,FIELD,DATA) /* ... */
#endif
//--------------------------------------------------
#ifdef TRACE_BLOCK
#   undef TRACE_BLOCK
#   define TRACE_BLOCK(MSG,BLOCK)                                       \
  CkPrintf ("%d TRACE_RESTART BLOCK       %s %s\n",CkMyPe(),BLOCK->name().c_str(), \
            std::string(MSG).c_str());                                  \
  fflush(stdout);
#   define TRACE_READER(MSG,READER)                                     \
  CkPrintf ("%d TRACE_RESTART READER    %s %d\n",                          \
            CkMyPe(),std::string(MSG).c_str(),READER->thisIndex);       \
  fflush(stdout);
#   define TRACE_SIMULATION(MSG,SIMULATION)                             \
  CkPrintf ("%d TRACE_RESTART SIMULATION %s %d\n",                      \
            CkMyPe(),std::string(MSG).c_str(),SIMULATION->thisIndex);   \
  fflush(stdout);
#else
#   define TRACE_BLOCK(MSG,BLOCK) /* ... */
#   define TRACE_READER(MSG,READER)  /* ... */
#   define TRACE_SIMULATION(MSG,READER)  /* ... */
#endif

//----------------------------------------------------------------------
// STARTUP
//----------------------------------------------------------------------

void Block::restart_enter_()
{
  TRACE_BLOCK("restart_enter_",this);
  const std::string restart_dir  = cello::config()->initial_restart_dir;
  if (index_.is_root()) {
    proxy_simulation[0].p_restart_enter(restart_dir);
  }
}

//----------------------------------------------------------------------

void Simulation::p_restart_enter (std::string name_dir)
{
  TRACE_SIMULATION("p_restart_enter_",this);
  // [ Called on root process only ]

  restart_directory_ = name_dir;

  // Open and read the checkpoint file_list file
  restart_stream_file_list_ = file_open_file_list_(restart_directory_);
  restart_stream_file_list_ >> restart_num_files_;

  // set synchronization
  TRACE_SYNC(sync_restart_created_,"sync_restart_created_ set_stop()");
  sync_restart_created_.set_stop(restart_num_files_);
  TRACE_SYNC(sync_restart_next_,"sync_restart_next_ set_stop()");
  sync_restart_next_.set_stop(restart_num_files_);

  // Create new empty IoEnzoReader chare array and distribute to other processing elements
  CProxy_MappingIo io_map  = CProxy_MappingIo::ckNew(restart_num_files_);

  CkArrayOptions opts(restart_num_files_);
  opts.setMap(io_map);
  // create array
  proxy_io_enzo_reader = CProxy_IoEnzoReader::ckNew(opts);
}

//----------------------------------------------------------------------

IoEnzoReader::IoEnzoReader()
  : CBase_IoEnzoReader(),
    name_dir_(),
    name_file_(),
    max_level_(),
    stream_block_list_(),
    file_(nullptr),
    sync_blocks_(),
    io_msg_check_(),
    level_(0),
    block_name_list_(),
    block_level_list_(),
    blocks_in_level_()
{
  
  proxy_enzo_simulation.p_io_reader_created();
}

//----------------------------------------------------------------------

void EnzoSimulation::p_io_reader_created()
{
  TRACE_SIMULATION("p_io_reader_created",this);
  // Wait for all io_readers to be created
  TRACE_SYNC(sync_restart_created_,"sync_restart_created_ next()");
  if (sync_restart_created_.next()) {
    // distribute array proxy to other simulation objects
    proxy_enzo_simulation.p_set_io_reader(proxy_io_enzo_reader);
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_set_io_reader(CProxy_IoEnzoReader io_enzo_reader)
{
  TRACE_SIMULATION("EnzoSimulation::p_set_io_reader()",this);
  proxy_io_enzo_reader = io_enzo_reader;
  CkCallback callback(CkIndex_Simulation::r_restart_start(NULL),0,
                      proxy_simulation);
  contribute(callback);
}

//----------------------------------------------------------------------

void Simulation::r_restart_start (CkReductionMsg * msg)
{
  TRACE_SIMULATION("EnzoSimulation::r_restart_start()",this);
  delete msg;
  // [ Called on root process only ]

  // Insert an IoEnzoReader element for each file
  const int max_level = cello::config()->mesh_max_level;
  for (int i=0; i<restart_num_files_; i++) {
    std::string restart_file;
    restart_stream_file_list_ >> restart_file;
    // Create ith io_reader to read name_file
    proxy_io_enzo_reader[i].p_init_root
      (restart_directory_,restart_file,max_level);
  }
}

//----------------------------------------------------------------------
// LEVEL 0
//----------------------------------------------------------------------

void IoEnzoReader::p_init_root
(std::string name_dir, std::string name_file, int max_level)
{
  TRACE_READER("p_init_root()",this);
  // save initialization parameters
  name_dir_  = name_dir;
  name_file_ = name_file;
  max_level_ = max_level;

  io_msg_check_.resize(max_level+1);

  stream_block_list_ = stream_open_blocks_(name_dir, name_file);

  // open the HDF5 file
  file_open_block_list_(name_dir,name_file);

  sync_blocks_.reset();
  TRACE_SYNC(sync_blocks_,"sync_blocks_ reset()");

  // Read global attributes
  file_read_hierarchy_();

  blocks_in_level_.resize(max_level+1);

  std::string block_name;
  int block_level;
  // Read list of blocks and associated refinement levels
  while (read_block_list_(block_name,block_level)) {

    // save block name and level
    block_name_list_.push_back(block_name);
    block_level_list_.push_back(block_level);

    if (block_level <= 0) {
      // count root-level blocks for synchronization
      // (including negative level blocks)
      ++ sync_blocks_;
      TRACE_SYNC(sync_blocks_,"sync_blocks_ inc_stop(1)");
    } else {
      // count blocks per refined level for array allocation below
      ++blocks_in_level_[block_level];
    }
  }

  // Allocate io_msg_check_ array
  io_msg_check_.resize(max_level+1);
  for (int level=1; level<=max_level; level++) {
    io_msg_check_[level].resize(blocks_in_level_[level]);
  }

  // initialize vector of array offsets into io_msg_check_ for each level
  std::vector<int> level_index;
  level_index.resize(max_level+1);
  std::fill(level_index.begin(), level_index.end(), 0);

  // Save block's meta-data and initialize root-level blocks
  for (int i=0; i<block_name_list_.size(); i++) {

    // Increment block counter so know when to alert pe=0 when done
    const int block_level = block_level_list_[i];
    std::string block_name = block_name_list_[i];
    int i_level = level_index[std::max(block_level,0)]++;

    // For each Block in the file, read the block data and create the
    // new Block

    EnzoMsgCheck * msg_check = new EnzoMsgCheck;

    // Save msg_check for later if block is in a refined level
    if (block_level > 0) {
      io_msg_check_[block_level][i_level] = msg_check;
    }

    // create new IoEnzoBlock object for root blocks, use io_msg_check_
    // element for refined blocks
    IoEnzoBlock * io_enzo_block = new IoEnzoBlock;

    file_read_block_ (msg_check, block_name, io_enzo_block);

    // save this file IoReader index
    msg_check->index_file_ = thisIndex;

    // get Block's index
    int v3[3];
    msg_check->io_block_->index(v3);
    Index index;
    index.set_values(v3);

    // if non-refined block, initialize
    if (block_level <= 0) {

      // Block exists--send its data
#ifdef DEBUG_RESTART
      msg_check->print("send");
      msg_check->data_msg_->print("send");
#endif
      enzo::block_array()[index].p_restart_set_data(msg_check);

    }
  }

  // close the HDF5 file
  file_close_block_list_();

  // self + 1
  ++ sync_blocks_;
  block_ready_();
}

//----------------------------------------------------------------------

void EnzoBlock::p_restart_set_data(EnzoMsgCheck * msg_check)
{
#ifdef DEBUG_RESTART
  msg_check->print("recv");
  msg_check->data_msg_->print("read");
#endif
  restart_set_data_ (msg_check);
}

void EnzoBlock::restart_set_data_(EnzoMsgCheck * msg_check)
{
  TRACE_BLOCK("EnzoBlock::restart_set_data()",this);
  const int index_file = msg_check->index_file_;
  msg_check->update(this);
  msg_check->get_adapt(adapt_);

  PRINT_FIELD_RESTART("recv","density",data());
  delete msg_check;
  proxy_io_enzo_reader[index_file].p_block_ready();
}

//----------------------------------------------------------------------

void IoEnzoReader::p_block_ready()
{ block_ready_(); }

void IoEnzoReader::block_ready_()
{
  TRACE_READER("[p_]block_ready()",this);
  TRACE_SYNC(sync_blocks_,"sync_blocks_ next()");
  // Wait for all of the reader's blocks in the current level to be
  // ready
  if (sync_blocks_.next()) {
    proxy_enzo_simulation[0].p_restart_next_level();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_restart_next_level()
{
  // [ Called on root process only ]
  TRACE_SIMULATION("EnzoSimulation::p_restart_next_level()",this);
  TRACE_SYNC(sync_restart_next_,"sync_restart_next_ next()");
  if (sync_restart_next_.next()) {
    const int max_level = cello::config()->mesh_max_level;
    if (++restart_level_ <= max_level) {
      proxy_io_enzo_reader.p_create_level(restart_level_);
    } else {
      enzo::block_array().doneInserting();
      enzo::block_array().p_restart_done();
    }
  }
}

//----------------------------------------------------------------------
// LEVEL K
//----------------------------------------------------------------------

void IoEnzoReader::p_create_level (int level)
{
  level_ = level;
  TRACE_READER("p_create_level()",this);
  const int num_blocks_level = io_msg_check_[level].size();
  sync_blocks_.reset();
  sync_blocks_.set_stop(num_blocks_level+1);
  for (int i=0; i<num_blocks_level; i++) {
    IoBlock * io_block = io_msg_check_[level][i]->io_block();
    int i3[3];
    io_block->index(i3);
    Index index;
    index.set_values(i3);
    Index index_parent = index.index_parent();
    int ic3[3];
    index.child(level,ic3,ic3+1,ic3+2);
    enzo::block_array()[index_parent].p_restart_refine(ic3,thisIndex);
  }
  // self
  block_created_();
}
//----------------------------------------------------------------------

void EnzoBlock::p_restart_refine(int ic3[3],int io_reader)
{
  TRACE_BLOCK("EnzoBlock::p_restart_refine()",this);
  FieldData * field_data = data()->field_data();

  int nx,ny,nz;
  field_data->size(&nx,&ny,&nz);

  Index index_child = index_.index_child(ic3);

  // Create FieldFace for interpolating field data to child ic3[]

  int if3[3] = {0,0,0};
  int g3[3];
  cello::field_descr()->ghost_depth(0,g3,g3+1,g3+2);
  Refresh * refresh = new Refresh;
  refresh->add_all_data();
  FieldFace * field_face = create_face
    (if3,ic3,g3, refresh_fine, refresh, true);

  // Create data message object to send
  DataMsg * data_msg = new DataMsg;

  bool is_new;
  data_msg -> set_field_face (field_face,is_new=false);
  data_msg -> set_field_data (field_data,is_new=false);

  const Factory * factory = cello::simulation()->factory();

  // Create the child object with interpolated data

  int narray = 0;
  char * array = 0;
  int num_field_data = 1;

  factory->create_block
    (
     data_msg,
     thisProxy, index_child,
     nx,ny,nz,
     num_field_data,
     adapt_step_,
     cycle_,time_,dt_,
     narray, array, refresh_fine,
     27,
     &child_face_level_curr_.data()[27*IC3(ic3)],
     &adapt_,
     cello::simulation(),
     io_reader);

  delete [] array;
  array = 0;

  children_.push_back(index_child);

  // }
  adapt_.set_valid(false);
  is_leaf_ = false;
  
}

//----------------------------------------------------------------------

void IoEnzoReader::p_block_created()
{ block_created_(); }

void IoEnzoReader::block_created_()
{
  TRACE_READER("p_block_created()",this);
  TRACE_SYNC(sync_blocks_,"sync_blocks_ next()");
  if (sync_blocks_.next()) {
    proxy_enzo_simulation[0].p_restart_level_created();
  }
}

//----------------------------------------------------------------------

void EnzoSimulation::p_restart_level_created()
{
  TRACE_SIMULATION("EnzoSimulation::p_restart_level_created()",this);
  TRACE_SYNC(sync_restart_created_,"sync_restart_created_ next()");
  if (sync_restart_created_.next()) {
    proxy_io_enzo_reader.p_init_level(restart_level_);
  }
}

//----------------------------------------------------------------------

void IoEnzoReader::p_init_level (int level)
{
  TRACE_READER("p_init_level()",this);
  const int num_blocks_level = io_msg_check_[level].size();
  TRACE_SYNC(sync_blocks_,"sync_blocks_ reset()");
  sync_blocks_.reset();
  TRACE_SYNC(sync_blocks_,"sync_blocks_ set_stop()");
  sync_blocks_.set_stop(num_blocks_level+1);
  // Loop through blocks in the given level
  for (int i=0; i<num_blocks_level; i++) {

    EnzoMsgCheck * msg_check = io_msg_check_[level][i];

    // Get the current Block's index
    IoBlock * io_block = msg_check->io_block();
    int i3[3];
    io_block->index(i3);
    Index index;
    index.set_values(i3);

    msg_check->index_file_ = thisIndex;
#ifdef DEBUG_RESTART
    msg_check->print("send");
    msg_check->data_msg_->print("send");
#endif
    enzo::block_array()[index].p_restart_set_data(msg_check);
  }
  // self + 1
  block_ready_();
}

//----------------------------------------------------------------------

void EnzoBlock::p_restart_done()
{
  TRACE_BLOCK("EnzoBlock::p_restart_done()",this);
  adapt_exit_();
}

//======================================================================

std::ifstream Simulation::file_open_file_list_(std::string name_dir)
{
  std::string name_file = name_dir + "/check.file_list";

  std::ifstream stream_file_list (name_file);

  ASSERT1("Simulation::file_copen_file_list_",
          "Cannot open hierarchy file %s for reading",
          name_file.c_str(),stream_file_list);

  return stream_file_list;
}

//----------------------------------------------------------------------

std::ifstream IoEnzoReader::stream_open_blocks_
(std::string name_dir, std::string name_file)
{
  std::string name_file_full = name_dir + "/" + name_file + ".block_list";

  std::ifstream stream_block_list (name_file_full);

  ASSERT1("Simulation::create_block_list_",
          "Cannot open block_list file %s for reading",
          name_file_full.c_str(),stream_block_list);

  return stream_block_list;
}

//----------------------------------------------------------------------

void IoEnzoReader::file_read_hierarchy_()
{
  // Simulation data
  IoSimulation io_simulation = (cello::simulation());
  for (size_t i=0; i<io_simulation.meta_count(); i++) {

    void * buffer;
    std::string name;
    int type_scalar;
    int nx,ny,nz;

    // Get object's ith metadata
    io_simulation.meta_value(i,& buffer, &name, &type_scalar, &nx,&ny,&nz);

    // Read object's ith metadata
    file_->file_read_meta(buffer,name.c_str(),&type_scalar,&nx,&ny,&nz);
  }

  io_simulation.save_to(cello::simulation());

  // Hierarchy data
  IoHierarchy io_hierarchy = (cello::hierarchy());
  for (size_t i=0; i<io_hierarchy.meta_count(); i++) {

    void * buffer;
    std::string name;
    int type_scalar;
    int nx,ny,nz;

    // Get object's ith metadata
    io_hierarchy.meta_value(i,& buffer, &name, &type_scalar, &nx,&ny,&nz);

    // Read object's ith metadata
    file_->file_read_meta(buffer,name.c_str(),&type_scalar,&nx,&ny,&nz);
  }
  io_hierarchy.save_to(cello::hierarchy());
}

//----------------------------------------------------------------------

void IoEnzoReader::file_read_block_
(EnzoMsgCheck * msg_check,
 std::string    name_block,
 IoEnzoBlock *  io_block)
{
  // Open HDF5 group for the block
  std::string group_name = "/" + name_block;
  file_->group_chdir(group_name);
  file_->group_open();

  // Read the Block's attributes
  msg_check->set_io_block(io_block);
  read_meta_(file_, io_block, "group");

  // Read block Adapt
  int type, size;
  file_->group_read_meta
    (msg_check->adapt_buffer_,"adapt_buffer",&type,&size);

  DataMsg * data_msg = new DataMsg;
  msg_check->data_msg_ = data_msg;

  // Create and allocate the data object
  int nx,ny,nz;
  int root_blocks[3];
  int root_size[3];
  double xm,ym,zm;
  double xp,yp,zp;
  Hierarchy * hierarchy = cello::hierarchy();
  hierarchy->root_blocks(root_blocks,root_blocks+1,root_blocks+2);
  hierarchy->root_size(root_size,root_size+1,root_size+2);
  hierarchy->lower(&xm,&ym,&zm);
  hierarchy->upper(&xp,&yp,&zp);
  nx=root_size[0]/root_blocks[0];
  ny=root_size[1]/root_blocks[1];
  nz=root_size[2]/root_blocks[2];

  int num_field_blocks = 1;
  FieldDescr * field_descr = cello::field_descr();
  ParticleDescr * particle_descr = cello::particle_descr();

  Data * data = new Data
    (nx, ny, nz, num_field_blocks, xm,xp, ym,yp, zm,zp,
     field_descr, particle_descr);

  data->allocate();

  // Loop through fields and read them in

  const int num_fields = field_descr->num_permanent();
  Field field = data->field();

  if (num_fields > 0) {
    // If any fields, add them to DataMsg
    Refresh * refresh = new Refresh;
    refresh->add_all_data();
    FieldFace  * field_face = new FieldFace(cello::rank());

    field_face -> set_refresh_type (refresh_same);
    field_face -> set_child (0,0,0);
    field_face -> set_face (0,0,0);
    field_face -> set_ghost(true,true,true);
    field_face -> set_refresh(refresh,true);
    bool is_new;
    data_msg -> set_field_face (field_face,        is_new=true);
    data_msg -> set_field_data (data->field_data(),is_new=true);
  }
  for (int i_f=0; i_f<num_fields; i_f++) {

    const std::string field_name = field_descr->field_name(i_f);
    int index_field = field_descr->field_id(field_name);

    const std::string dataset_name = std::string("field_") + field_name;
    int m4[4];
    int type_data = type_unknown;
    file_->data_open (dataset_name, &type_data,
                      m4,m4+1,m4+2,m4+3);
    int mx,my,mz;
    int gx,gy,gz;

    field.dimensions(index_field,&mx,&my,&mz);
    field.ghost_depth(index_field,&gx,&gy,&gz);

    double lower[3];
    double upper[3];
    io_block->lower(lower);
    io_block->upper(upper);
    char * buffer = field.values(field_name);

    file_read_dataset_
      (buffer, type_data, mx,my,mz,m4);

    file_->data_close();

  }

  // Read in particle data

  Particle particle = data->particle();

  // for each particle type
  const int num_types = particle_descr->num_types();
  if (num_types > 0) {
    // If any fields, add them to DataMsg
    data_msg -> set_particle_data (particle.particle_data(),true);
  }
  for (int it=0; it<num_types; it++) {

    const std::string particle_name = particle_descr->type_name(it);
    // for each attribute of the particle type
    const int na = particle_descr->num_attributes(it);
    for (int ia=0; ia<na; ia++) {
      const std::string attribute_name = particle_descr->attribute_name(it,ia);

      const std::string dataset_name =
        std::string("particle_") + particle_name + "_" + attribute_name;
      int m4[4];
      int type_data = type_unknown;
      file_->data_open (dataset_name, &type_data,
                        m4,m4+1,m4+2,m4+3);

      const int np = m4[0];

      // allocate particles at start once count is known
#ifdef DEBUG_RESTART
      CkPrintf ("DEBUG_RESTART num_particles %s %d:%d %d\n",
                name_block.c_str(),it,ia,np);
#endif
      if (ia==0) {
        particle.insert_particles(it,np);
      }

      // read particle attribute into single array first...
      union {
        char * buffer;
        float * buffer_single;
        double * buffer_double;
        int8_t * buffer_int8;
        int16_t * buffer_int16;
        int32_t * buffer_int32;
        int64_t * buffer_int64;
      };

      buffer = file_->allocate_buffer(np,type_data);

      int nx=m4[0];
      int ny=m4[1];
      int nz=m4[2];
      file_read_dataset_(buffer, type_data, nx,ny,nz,m4);

      // ...then copy to particle batches

      if (type_data == type_single) {
        copy_buffer_to_particle_attribute_
          (buffer_single, particle, it, ia, np);
      } else if (type_data == type_double) {
        copy_buffer_to_particle_attribute_
          (buffer_double, particle, it, ia, np);
      } else if (type_data == type_int8) {
        copy_buffer_to_particle_attribute_
          (buffer_int8, particle, it, ia, np);
      } else if (type_data == type_int16) {
        copy_buffer_to_particle_attribute_
          (buffer_int16, particle, it, ia, np);
      } else if (type_data == type_int32) {
        copy_buffer_to_particle_attribute_
          (buffer_int32, particle, it, ia, np);
      } else if (type_data == type_int64) {
        copy_buffer_to_particle_attribute_
          (buffer_int64, particle, it, ia, np);
      } else {
        ERROR1 ("IoEnzoReader::file_read_block_()",
                "Unsupported particle type_data %d",
                type_data);
      }
      delete [] buffer;
    }
  }
  file_->group_close();
}

//----------------------------------------------------------------------

template <class T>
void IoEnzoReader::copy_buffer_to_particle_attribute_
(T * buffer, Particle particle, int it, int ia, int np)
{
  for (int ip=0; ip<np; ip++) {
    int ib,io;
    particle.index(ip,&ib,&io);
    T * batch = (T *) particle.attribute_array(it,ia,ib);
    batch[io] = buffer[ip];
  }
}

//----------------------------------------------------------------------

void IoEnzoReader::read_meta_
( FileHdf5 * file, Io * io, std::string type_meta )
{
  for (size_t i=0; i<io->meta_count(); i++) {

    void * buffer;
    std::string name;
    int type_scalar;
    int nx,ny,nz;

    // Get object's ith metadata
    io->meta_value(i,& buffer, &name, &type_scalar, &nx,&ny,&nz);

    // Read object's ith metadata
    if ( type_meta == "group" ) {
      file->group_read_meta(buffer,name.c_str(),&type_scalar,&nx,&ny,&nz);
    } else if (type_meta == "file") {
      file->file_read_meta(buffer,name.c_str(),&type_scalar,&nx,&ny,&nz);
    } else {
      ERROR1 ("MethodOutput::read_meta_()",
              "Unknown type_meta \"%s\"",
              type_meta.c_str());
    }
  }
}
//----------------------------------------------------------------------

bool IoEnzoReader::read_block_list_(std::string & block_name, int & level)
{
  bool value (stream_block_list_ >> block_name >> level);
  return value;
}

//----------------------------------------------------------------------
void IoEnzoReader::file_open_block_list_
(std::string path_name, std::string file_name)
{
  // Create File
  file_name = file_name + ".h5";
  file_ = new FileHdf5 (path_name, file_name);
  file_->file_open();
}

//----------------------------------------------------------------------

void IoEnzoReader::file_read_dataset_
(char * buffer, int type_data,
 int nx, int ny, int nz,
 int m4[4])
{
  // Read the domain dimensions

  // field size
  int n4[4];
  n4[0] = n4[1] = n4[2] = n4[3] = 1;
  n4[0] = nx;
  n4[1] = ny;
  n4[2] = nz;

  // determine offsets
  int o4[4] = {0,0,0,0};

  // open the dataspace
  file_-> data_slice
    (m4[0],m4[1],m4[2],m4[3],
     n4[0],n4[1],n4[2],n4[3],
     o4[0],o4[1],o4[2],o4[3]);

  // create memory space
  file_->mem_create (nx,ny,nz,nx,ny,nz,0,0,0);

  file_->data_read (buffer);
}
//----------------------------------------------------------------------

void IoEnzoReader::file_close_block_list_()
{
  file_->data_close();
  file_->file_close();
  delete file_;
}


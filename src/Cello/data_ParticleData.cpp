// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Aug 15 11:48:15 PDT 2014
/// @brief    Implementation of the ParticleData class

#include "data.hpp"

//----------------------------------------------------------------------

ParticleData::ParticleData()
  : attribute_array_(),
    attribute_align_(),
    particle_count_()
{
}

//----------------------------------------------------------------------

void ParticleData::pup (PUP::er &p)
{
  p | attribute_array_;
  p | attribute_align_;
  p | particle_count_;
}

//----------------------------------------------------------------------

char * ParticleData::attribute_array (ParticleDescr * particle_descr,
				      int it,int ib,int ia)
{

  bool in_range = true;
  if ( !(0 <= it && it < particle_descr->num_types()) )
    in_range = false;
  if ( !(0 <= ib && ib < num_batches(it)) )
    in_range = false;
  if ( !(0 <= ia && ia < particle_descr->num_attributes(it)) )
    in_range = false;

  char * array = NULL;
  if ( in_range ) {
    int offset = particle_descr->attribute_offset(it,ia);
    int align =  attribute_align_[it][ib];
    array = &attribute_array_[it][ib][0] + (offset + align);
  }
  return array;
}

//----------------------------------------------------------------------

int ParticleData::num_batches (int it) const
{
  if ( !(0 <= it && it < (int) attribute_array_.size() )) return 0;

  return attribute_array_[it].size();
}

//----------------------------------------------------------------------

int ParticleData::num_particles 
(ParticleDescr * particle_descr, int it, int ib) const
{
  if ( !(0 <= it && it < particle_descr->num_types()) ) return 0;
  if ( !(0 <= ib && ib < num_batches(it)) ) return 0;

  return particle_count_[it][ib];
}

//----------------------------------------------------------------------

int ParticleData::num_particles 
(ParticleDescr * particle_descr,int it) const
{
  int nb = num_batches(it);
  int np = 0;
  for (int ib=0; ib<nb; ib++) {
    np += num_particles(particle_descr,it,ib);
  }
  return np;
}

//----------------------------------------------------------------------

int ParticleData::num_particles 
(ParticleDescr * particle_descr) const
{
  const int nt = particle_descr->num_types();
  int np = 0;
  for (int it=0; it<nt; it++) {
    np += num_particles(particle_descr,it);
  }
  return np;
}

//----------------------------------------------------------------------

int ParticleData::insert_particles 
(ParticleDescr * particle_descr,
 int it, int np)
{

  check_arrays_(particle_descr,__FILE__,__LINE__);

  // find indices of last batch and particle in batch for return value
  int ib_last = num_batches(it) - 1;
  int ip_last = num_particles(particle_descr,it,ib_last);

  int np_left = np;

  const int np_batch = particle_descr->batch_size();

  if (ib_last < 0 || ip_last == np_batch) {
    ib_last++;
    ip_last = 0;
  }
  int ib_this = ib_last;
  int ip_start = ip_last;
  while (np_left > 0) {

    // number of particles to add in this batch
    int np_this = std::min(np_batch,np_left) - ip_start;

    // resize arrays for a new batch if needed
    if ( ib_this > num_batches(it) - 1) {
      attribute_array_[it].resize(ib_this+1);
      attribute_align_[it].resize(ib_this+1);
      particle_count_[it].resize(ib_this+1);
    }

    // allocate particles
    
    resize_array_(particle_descr,it,ib_this,ip_start+np_this);

    // prepare for next batch

    np_left -= np_this;
    ib_this++;
    ip_start=0;
  }

  // return global index of first particle
  return (ip_last + np_batch*ib_last);
}

//----------------------------------------------------------------------

void ParticleData::delete_particles 
(ParticleDescr * particle_descr,
 int it, int ib, const bool * mask)
{
  check_arrays_(particle_descr,__FILE__,__LINE__);

  const bool interleaved = particle_descr->interleaved(it);

  int npd=0;

  const int na = particle_descr->num_attributes(it);
  const int np = num_particles(particle_descr,it,ib);
  int mp = particle_descr->particle_bytes(it);

  for (int ip=0; ip<np; ip++) {
    if (mask[ip]) {
      npd++;
    } else if (npd>0) {
      for (int ia=0; ia<na; ia++) {
	if (!interleaved) 
	  mp = particle_descr->attribute_bytes(it,ia);
	int ny = particle_descr->attribute_bytes(it,ia);
	char * a = attribute_array(particle_descr,it,ib,ia);
	for (int iy=0; iy<ny; iy++) {
	  a [iy + mp*(ip-npd)] = a [iy + mp*ip];
	}
      }
    }
  }

  if (npd>0) {
    resize_array_(particle_descr,it,ib,np-npd);
  }
}

//----------------------------------------------------------------------

void ParticleData::split_particles 
(ParticleDescr * particle_descr, 
 int it, int ib, const bool *mask,
 ParticleData * particle_data_dest)
{
  check_arrays_(particle_descr,__FILE__,__LINE__);
  const int na = particle_descr->num_attributes(it);
  const int np = num_particles(particle_descr,it,ib);

  // Return if no particles deleted
  int nd=0;
  for (int ip=0; ip<np; ip++) {
    if (mask[ip]) nd++;
  }
  if (nd==0) return;

  // allocate nd particles

  int j = particle_data_dest->insert_particles(particle_descr,it,nd);
  int jb,jp;
  particle_descr->index(j,&jb,&jp);

  // copy particles to be deleted

  const bool interleaved = particle_descr->interleaved(it);
  const int mb = particle_descr->batch_size();

  int mp = particle_descr->particle_bytes(it);

  for (int ip=0; ip<np; ip++) {
    if (mask[ip]) {
      for (int ia=0; ia<na; ia++) {
	if (!interleaved) 
	  mp = particle_descr->attribute_bytes(it,ia);
	int ny = particle_descr->attribute_bytes(it,ia);
	char * a_src = attribute_array(particle_descr,it,ib,ia);
	char * a_dst = attribute_array(particle_descr,it,jb,ia);
	for (int iy=0; iy<ny; iy++) {
	  a_dst [iy + mp*jp] = a_src [iy + mp*ip];
	}
      }
      jp = (jp+1) % mb;
      if (jp==0) jb++;
    }
  }

  // delete particles
  delete_particles(particle_descr,it,ib,mask);
}

//----------------------------------------------------------------------

void ParticleData::compress (ParticleDescr * particle_descr)
{
  const int nt = particle_descr->num_types();
  for (int it=0; it<nt; it++) {
    compress (particle_descr,it);
  }
}

//----------------------------------------------------------------------

void ParticleData::compress (ParticleDescr * particle_descr, int it)
{
  const int nb = num_batches(it);
  const int mb = particle_descr->batch_size();
  const int na = particle_descr->num_attributes(it);

  const bool interleaved = particle_descr->interleaved(it);

  // find first batch with space in it

  int ibs, ips; // source batch and particle indices
  int ibd, ipd; // destination batch and particle indices

  int npd; // number of particles in ibd batch
  int nps; // number of particles in ibs batch

  // find destination: first empty spot 

  ibd=0;
  npd=(ibd<nb) ? num_particles(particle_descr,it,ibd) : 0;
  while (ibd<nb && npd == mb) {
    ibd++;
  } // assert ibd == nb || npd < mb
  ipd = npd;

  // first source: next non-empty spot
  ibs = ibd + 1;
  ips = 0;
  nps = num_particles(particle_descr,it,ibs);

  resize_array_ (particle_descr,it,ibd,mb);
  npd = mb;

  int mp = particle_descr->particle_bytes(it);

  while (ibs < nb && ips < nps) {

    for (int ia=0; ia<na; ia++) {
      if (!interleaved) {
	mp = particle_descr->attribute_bytes(it,ia);
      }
      int ny = particle_descr->attribute_bytes(it,ia);
      char * as = attribute_array(particle_descr,it,ibs,ia);
      char * ad = attribute_array(particle_descr,it,ibd,ia);
      for (int iy=0; iy<ny; iy++) {
	ad [iy + mp*ipd] = as [iy + mp*ips];
      }
    }

    ipd++;
    if (ipd>=npd) {
      ipd = 0;
      ibd++;
      npd = mb;
      if (ibd < nb) resize_array_ (particle_descr,it,ibd,mb);
    }

    ips++;
    if (ips>=nps) {
      ips = 0;
      ibs++;
      if (ibs < nb) nps = num_particles(particle_descr,it,ibs);
    }

  }
}
  

//----------------------------------------------------------------------

float ParticleData::efficiency (ParticleDescr * particle_descr)
{
  long bytes_min=0,bytes_used=0;
  const int mb = particle_descr->batch_size();

  const int nt = particle_descr->num_types();
  for (int it=0; it<nt; it++) {
    const int nb = num_batches(it);
    const int mp = particle_descr->particle_bytes(it);
    for (int ib=0; ib<nb; ib++) {
      const int np = num_particles(particle_descr,it,ib);
      bytes_min += np*mp;
      bytes_used += mb*mp;
    }
  }
  return 1.0*bytes_min/bytes_used;

}

//----------------------------------------------------------------------

float ParticleData::efficiency (ParticleDescr * particle_descr, int it)
{
  long bytes_min=0,bytes_used=0;
  const int mb = particle_descr->batch_size();

  const int nb = num_batches(it);
  const int mp = particle_descr->particle_bytes(it);
  for (int ib=0; ib<nb; ib++) {
    const int np = num_particles(particle_descr,it,ib);
    bytes_min += np*mp;
    bytes_used += mb*mp;
  }
  return 1.0*bytes_min/bytes_used;
}

//----------------------------------------------------------------------

float ParticleData::efficiency (ParticleDescr * particle_descr, int it, int ib)
{

  const int mp = particle_descr->particle_bytes(it);
  const int np = num_particles(particle_descr,it,ib);
  const int mb = particle_descr->batch_size();

  const long bytes_min  = np*mp;
  const long bytes_used = mb*mp;

  return 1.0*bytes_min/bytes_used;
  
}

//======================================================================


void ParticleData::resize_array_(ParticleDescr * particle_descr,
				 int it, int ib, int np)
{
  // store number of particles allocated
  particle_count_[it][ib] = np;

  const int mp = particle_descr->particle_bytes(it);

  if (!particle_descr->interleaved(it)) {
    np = particle_descr->batch_size();
  }
  attribute_array_[it][ib].resize(mp*(np) + (PARTICLE_ALIGN - 1) );
  char * array = &attribute_array_[it][ib][0];
  uintptr_t iarray = (uintptr_t) array;
  int defect = (iarray % PARTICLE_ALIGN);
  attribute_align_[it][ib] = (defect == 0) ? 0 : PARTICLE_ALIGN-defect;
}

//----------------------------------------------------------------------

void ParticleData::check_arrays_ (ParticleDescr * particle_descr,
		    std::string file, int line) const
{
  size_t nt = particle_descr->num_types();
  ASSERT4 ("ParticleData::check_arrays_",
	  "%s:%d attribute_array_ is size %d < %d",
	   file.c_str(),line,
	   attribute_array_.size(),nt,
	   attribute_array_.size()>=nt);
  ASSERT4 ("ParticleData::check_arrays_",
	  "%s:%d particle_count_ is size %d < %d",
	   file.c_str(),line,
	   particle_count_.size(),nt,
	   particle_count_.size()>=nt);
  ASSERT4 ("ParticleData::check_arrays_",
	  "%s:%d attribute_align_ is size %d < %d",
	   file.c_str(),line,
	   attribute_align_.size(),nt,
	   attribute_align_.size()>=nt);

  for (size_t it=0; it<nt; it++) {
    size_t nb = num_batches(it);

    ASSERT5 ("ParticleData::check_arrays_",
	     "%s:%d attribute_array_[%d] is size %d < %d",
	     file.c_str(),line,
	     it,attribute_array_[it].size(),nb,
	     attribute_array_[it].size()>=nb);
    ASSERT5 ("ParticleData::check_arrays_",
	     "%s:%d particle_count_[%d] is size %d < %d",
	     file.c_str(),line,
	     it,particle_count_[it].size(),nb,
	     particle_count_[it].size()>=nb);
    ASSERT5 ("ParticleData::check_arrays_",
	     "%s:%d attribute_align_[%d] is size %d < %d",
	     file.c_str(),line,
	     it,attribute_align_[it].size(),nb,
	     attribute_align_[it].size()>=nb);
  }
}
// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "mumps.h"
#include "../trace.h"
#include "../error.h"
#include "../utils.h"
#include "../callstack.h"

#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  #define MUMPS         dmumps_c
  #define MUMPS_STRUCT  DMUMPS_STRUC_C
#else
  #define MUMPS         zmumps_c
  #define MUMPS_STRUCT  ZMUMPS_STRUC_C
#endif

#define USE_COMM_WORLD  -987654

#ifdef WITH_MUMPS

extern "C" {
  extern void MUMPS(MUMPS_STRUCT *idptr);
}

#else
#endif

// Binary search for the location of a particular CSC/CSR matrix entry.
//
// Typically, we search for the index into Ax that corresponds to a given 
// row (CSC) or column (CSR) ('idx') among indices of nonzero values in 
// a particular column (CSC) or row (CSR) ('Ai').
//
static int find_position(int *Ai, int Alen, int idx) {
  _F_
  assert (idx >= 0);
  
  register int lo = 0, hi = Alen - 1, mid;

  while (1) 
  {
    mid = (lo + hi) >> 1;

    if (idx < Ai[mid]) hi = mid - 1;
    else if (idx > Ai[mid]) lo = mid + 1;
    else break;

    // Sparse matrix entry not found (raise an error when trying to add 
    // value to this position, return 0 when obtaining value there).
    if (lo > hi) mid = -1;
  }
  return mid;
}

MumpsMatrix::MumpsMatrix()
{
  _F_
  nnz = 0;
  size = 0;
  irn = NULL;
  jcn = NULL;
  Ax = NULL;
  Ap = NULL;
  Ai = NULL;
}

MumpsMatrix::~MumpsMatrix()
{
  _F_
  free();
}

void MumpsMatrix::alloc()
{
  _F_
  assert(pages != NULL);
  assert(size > 0);

  // initialize the arrays Ap and Ai
  Ap = new int [size + 1];
  MEM_CHECK(Ap);
  int aisize = get_num_indices();
  Ai = new int [aisize];
  MEM_CHECK(Ai);

  // sort the indices and remove duplicities, insert into Ai
  int i, pos = 0;
  for (i = 0; i < size; i++) {
    Ap[i] = pos;
    pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
  }
  Ap[i] = pos;

  delete [] pages;
  pages = NULL;

  nnz = Ap[size];

  Ax = new mumps_scalar[nnz];
  memset(Ax, 0, sizeof(mumps_scalar) * nnz);

  irn = new int[nnz];
  memset(irn, 0, sizeof(int) * nnz);
  jcn = new int[nnz];
  memset(jcn, 0, sizeof(int) * nnz);
}

void MumpsMatrix::free()
{
  _F_
  nnz = 0;
  delete[] Ap; Ap = NULL;
  delete[] Ai; Ai = NULL;
  delete[] Ax; Ax = NULL;
  delete[] irn; irn = NULL;
  delete[] jcn; jcn = NULL;
}

scalar MumpsMatrix::get(int m, int n)
{
  _F_
  // Find m-th row in the n-th column.
  int mid = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
  // Return 0 if the entry has not been found.
  if (mid < 0) return 0.0;
  // Otherwise, add offset to the n-th column and return the value.
  if (mid >= 0) mid += Ap[n];
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  return Ax[mid];
#else
  return cplx(Ax[mid].r, Ax[mid].i);
#endif
}

void MumpsMatrix::zero()
{
  _F_
  memset(Ax, 0, sizeof(mumps_scalar) * Ap[size]);
}

void MumpsMatrix::add(int m, int n, scalar v)
{
  _F_
  // WARNING: The additional condition v != 0.0 used in (Pardiso/Umfpack)Matrix
  //          produced an error in neutronics-2-group-adapt (although tutorial-07
  //          ran well).
  if (m >= 0 && n >= 0) // ignore dirichlet DOFs
  {   
    // Find m-th row in the n-th column.
    int pos = find_position(Ai + Ap[n], Ap[n + 1] - Ap[n], m);
    // Make sure we are adding to an existing non-zero entry.
    if (pos < 0) 
      error("Sparse matrix entry not found");
    // Add offset to the n-th column.
    pos += Ap[n];
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    Ax[pos] += v;
#else
    Ax[pos].r += v.real();
    Ax[pos].i += v.imag();
#endif
    irn[pos] = m + 1;  // MUMPS is indexing from 1
    jcn[pos] = n + 1;
  }
}

void MumpsMatrix::add(int m, int n, scalar **mat, int *rows, int *cols)
{
  _F_
  for (int i = 0; i < m; i++)       // rows
    for (int j = 0; j < n; j++)     // cols
      add(rows[i], cols[j], mat[i][j]);
}

/// dumping matrix and right-hand side
///
bool MumpsMatrix::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  // TODO
  switch (fmt) 
  {
    case DF_NATIVE:
    case DF_PLAIN_ASCII:
      fprintf(file, "%d\n", size);
      fprintf(file, "%d\n", nnz);
      for (int i = 0; i < nnz; i++)
        fprintf(file, "%d %d " SCALAR_FMT "\n", irn[i], jcn[i], MUMPS_SCALAR(Ax[i]));
      return true;

    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", size, size, Ap[size], Ap[size]);
      for (int j = 0; j < size; j++)
        for (int i = Ap[j]; i < Ap[j + 1]; i++)
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)          
          fprintf(file, "%d %d " SCALAR_FMT "\n", Ai[i] + 1, j + 1, MUMPS_SCALAR(Ax[i]));
#else          
          fprintf(file, "%d %d %lf+%lfi\n", Ai[i] + 1, j + 1, MUMPS_SCALAR(Ax[i]));
#endif          
      fprintf(file, "];\n%s = spconvert(temp);\n", var_name);

      return true;

    case DF_HERMES_BIN: 
    {
      hermes_fwrite("H3DX\001\000\000\000", 1, 8, file);
      int ssize = sizeof(scalar);
      hermes_fwrite(&ssize, sizeof(int), 1, file);
      hermes_fwrite(&size, sizeof(int), 1, file);
      hermes_fwrite(&nnz, sizeof(int), 1, file);
      hermes_fwrite(Ap, sizeof(int), size + 1, file);
      hermes_fwrite(Ai, sizeof(int), nnz, file);
      hermes_fwrite(Ax, sizeof(mumps_scalar), nnz, file);
      return true;
    }

    default:
      return false;
  }
}

int MumpsMatrix::get_matrix_size() const
{
  _F_
  /*           Ax               Ai                 Ap                 */
  return (sizeof(scalar) + sizeof(int)) * nnz + sizeof(int)*(size+1)
          + 2 * sizeof(int) * nnz + sizeof(int);
  /*          irn, jcn                  nnz                           */    
}

double MumpsMatrix::get_fill_in() const
{
  _F_
  return Ap[size] / (double) (size * size);
}

// MumpsVector /////////////////////////////////////////////////////////////////////////////////////

MumpsVector::MumpsVector()
{
  _F_
  v = NULL;
  size = 0;
}

MumpsVector::~MumpsVector()
{
  _F_
  free();
}

void MumpsVector::alloc(int n)
{
  _F_
  free();
  size = n;
  v = new mumps_scalar[n];
  zero();
}

void MumpsVector::zero()
{
  _F_
  memset(v, 0, size * sizeof(mumps_scalar));
}

void MumpsVector::free()
{
  _F_
  delete [] v;
  v = NULL;
  size = 0;
}

void MumpsVector::set(int idx, scalar y)
{
  _F_
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  if (idx >= 0) v[idx] = y;
#else
  if (idx >= 0) {
    v[idx].r = y.real();
    v[idx].i = y.imag();
  }
#endif
}

void MumpsVector::add(int idx, scalar y)
{
  _F_
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
  if (idx >= 0) v[idx] += y;
#else
  if (idx >= 0) {
    v[idx].r += y.real();
    v[idx].i += y.imag();
  }
#endif
}

void MumpsVector::add(int n, int *idx, scalar *y)
{
  _F_
  for (int i = 0; i < n; i++)
    if (idx[i] >= 0) {
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
      v[idx[i]] += y[i];
#else
      v[idx[i]].r += y[i].real();
      v[idx[i]].i += y[i].imag();
#endif
    }
}

bool MumpsVector::dump(FILE *file, const char *var_name, EMatrixDumpFormat fmt)
{
  _F_
  switch (fmt) 
  {
    case DF_NATIVE:
    case DF_PLAIN_ASCII:
      for (int i = 0; i < size; i++)
        fprintf(file, SCALAR_FMT "\n", MUMPS_SCALAR(v[i]));

      return true;

    case DF_MATLAB_SPARSE:
      fprintf(file, "%% Size: %dx1\n%s = [\n", size, var_name);
      for (int i = 0; i < size; i++)
        fprintf(file, SCALAR_FMT "\n", MUMPS_SCALAR(v[i]));
      fprintf(file, " ];\n");
      return true;

    case DF_HERMES_BIN: 
    {
      hermes_fwrite("H3DR\001\000\000\000", 1, 8, file);
      int ssize = sizeof(scalar);
      hermes_fwrite(&ssize, sizeof(int), 1, file);
      hermes_fwrite(&size, sizeof(int), 1, file);
      hermes_fwrite(v, sizeof(scalar), size, file);
      return true;
    }

    default:
      return false;
  }
}

// MUMPS solver ////////////////////////////////////////////////////////////////////////////////////

MumpsSolver::MumpsSolver(MumpsMatrix *m, MumpsVector *rhs) :
  LinearSolver(), m(m), rhs(rhs)
{
  _F_
#ifdef WITH_MUMPS
#else
  error(MUMPS_NOT_COMPILED);
#endif
}

MumpsSolver::~MumpsSolver()
{
  _F_
#ifdef WITH_MUMPS
  //if (m != NULL) delete m;
  //if (rhs != NULL) delete rhs;
#endif
}

#ifdef WITH_MUMPS

// macro s.t. indices match Fortran documentation
#define ICNTL(I)            icntl[(I)-1]
#define MUMPS_INFO(id, I)   id->infog[(I)-1]
#define INFOG(I)            infog[(I)-1]

#define JOB_INIT            -1
#define JOB_END             -2

static bool check_status(MUMPS_STRUCT *id)
{
  _F_
  switch (id->INFOG(1)) {
    case 0: return true; // no error
    case -1: warning("Error occured on processor %d", MUMPS_INFO(id, 2)); break;
    // TODO: add the rest according to the MUMPS docs
    default: warning("INFOG(1) = %d", id->INFOG(1)); break;
  }
  return false;
}

#endif

bool MumpsSolver::solve()
{
  _F_
#ifdef WITH_MUMPS
  bool ret = false;
  assert(m != NULL);
  assert(rhs != NULL);

  TimePeriod tmr;

  MUMPS_STRUCT id;

  // Initialize a MUMPS instance
  id.job = JOB_INIT;
  id.par = 1; // host also performs calculations
  id.sym = 0; // 0 = unsymmetric
  id.comm_fortran=USE_COMM_WORLD;
  
  MUMPS(&id);
  check_status(&id);

  // matrix
  id.n = m->size;
  id.nz = m->nnz;
  id.irn = m->irn;
  id.jcn = m->jcn;
  id.a = m->Ax;

  // right-hand side
  id.rhs = new mumps_scalar[m->size];
  memcpy(id.rhs, rhs->v, m->size * sizeof(mumps_scalar));

  // No printings
  id.ICNTL(1) = -1;
  id.ICNTL(2) = -1;
  id.ICNTL(3) = -1;
  id.ICNTL(4) = 0;

  id.ICNTL(20) = 0; // centralized dense RHS
  id.ICNTL(21) = 0; // centralized dense solution
  
  id.job = 6; // 6 means Analysis + factorization + solve (FIXME: remove magic constant - see MUMPS docs)
  MUMPS(&id);
  
  ret = check_status(&id);

  if (ret) {
    delete [] sln;
    sln = new scalar[m->size];
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    for (int i = 0; i < rhs->size; i++)
      sln[i] = id.rhs[i];
#else
    for (int i = 0; i < rhs->size; i++)
      sln[i] = cplx(id.rhs[i].r, id.rhs[i].i);
#endif
  }

  // Terminate/free current instance
  id.job = JOB_END;
  MUMPS(&id);

  tmr.tick();
  time = tmr.accumulated();

  delete [] id.rhs;

  return ret;
#else
  return false;
#endif
}

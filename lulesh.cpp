#include "lulesh.h"
#include <sys/time.h>

double gettime() {
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return double(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

/******************************************************************************
 * tag distribution
 */
int iteration_tag_gen_step::execute(const int & t, lulesh_context & c) const {
  const int iteration = t;
  int num_blocks = c.num_blocks;
  int i,j,k;
  // Consume
  double previous_delta_time, delta_time, elapsed_time, targetdt, clock_time;
  c.delta_time.get(iteration, previous_delta_time);
  c.elapsed_time.get(iteration, elapsed_time);
  c.clock_time.get(iteration, clock_time);

  // Constants
  const double max_delta_time = c.max_delta_time;
  const double deltatimemultlb = 1.1;
  const double deltatimemultub = 1.2;
  double gnewdt = 1.0e+20;
  targetdt = c.stop_time - elapsed_time;

  if(iteration != 0) {
    double min_courant = 1e20;
    double min_hydro = 1e20;
    for (i=0; i < num_blocks; ++i) {
      for (j=0; j < num_blocks; ++j) {
        for (k=0; k < num_blocks; ++k) {
          double courant, hydro;
          c.mincourant.get(Quad(iteration, i, j, k), courant);
          c.minhydro.get(Quad(iteration, i, j, k), hydro);
          min_courant = min(courant, min_courant);
          min_hydro = min(hydro, min_hydro);
        }
      }
    }

    if (min_courant < gnewdt)
      gnewdt = min_courant / 2.0;

    if (min_hydro < gnewdt)
      gnewdt = min_hydro * 2.0 / 3.0;

    delta_time = gnewdt;
    double ratio = delta_time / previous_delta_time;

    if (ratio >= 1.0) {
      if (ratio < deltatimemultlb) {
        delta_time = previous_delta_time ;
      } else if (ratio > deltatimemultub) {
        delta_time= previous_delta_time * deltatimemultub;
      }
    }

    if (delta_time > max_delta_time)
      delta_time = max_delta_time;

  } else {
    delta_time = previous_delta_time;
  }

  // TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE
  if ((targetdt > delta_time) && (targetdt < (4.0 * delta_time / 3.0)) )
    targetdt = 2.0 * delta_time / 3.0;

  if (targetdt < delta_time)
    delta_time = targetdt ;

  // std::cout << "olddt = " << previous_delta_time << ", newdt = " << delta_time << "\n";

  if(elapsed_time < c.stop_time && (iteration < c.max_iter)) {

    elapsed_time += delta_time;
    double end_time = gettime();
    // Produce
    c.delta_time.put(iteration + 1, delta_time);
    c.elapsed_time.put(iteration + 1, elapsed_time);
    c.clock_time.put(iteration + 1, end_time);
    // std::cout << "iteration " << iteration << ", took: " << end_time - clock_time << "\n";
    if (c.debug) std::cout <<  end_time - clock_time << "\n";

    for (i=0; i < num_blocks; ++i) {
      for (j=0; j < num_blocks; ++j) {
        for (k=0; k < num_blocks; ++k) {
          c.node_tags.put(Quad(iteration+1, i, j, k));
          c.element_tags.put(Quad(iteration+1, i, j, k));
        }
      }
    }
    c.iteration.put(iteration + 1);
  }
  return CnC::CNC_Success;
}

// Computes Node Data for each: Element ID based tags
int compute_element_forces::execute(const Quad &q, lulesh_context & c) const {
  // Element ID based tags
  const int iteration = q[0];
  const int block_dim = c.block_dim;
  const int d2 = block_dim+1; // need extra space to cover node neighbors
  int xlow, ylow, zlow;
  int i,j,k;

  xlow = block_dim*q[1];
  ylow = block_dim*q[2];
  zlow = block_dim*q[3];

  int element_id, node_id, local_node_id, local_element_id;

  std::shared_ptr<Tile3d<double> > pressure_blk, viscosity_blk, volume_blk, ss_blk;
  std::shared_ptr<Tile3d<double> > posx, posy, posz, velx, vely, velz;

  // Get element data
  c.pressure.get(Quad(q[0]-1, q[1], q[2], q[3]), pressure_blk);
  c.viscosity.get(Quad(q[0]-1, q[1], q[2], q[3]), viscosity_blk);
  c.volume.get(Quad(q[0]-1, q[1], q[2], q[3]), volume_blk);
  c.sound_speed.get(Quad(q[0]-1, q[1], q[2], q[3]), ss_blk);

  // Get my node tile
  c.posx.get(Quad(q[0]-1, q[1], q[2], q[3]), posx);
  c.posy.get(Quad(q[0]-1, q[1], q[2], q[3]), posy);
  c.posz.get(Quad(q[0]-1, q[1], q[2], q[3]), posz);
  c.velx.get(Quad(q[0]-1, q[1], q[2], q[3]), velx);
  c.vely.get(Quad(q[0]-1, q[1], q[2], q[3]), vely);
  c.velz.get(Quad(q[0]-1, q[1], q[2], q[3]), velz);

  std::shared_ptr<Tile3d<double> > forcesx =
    std::make_shared<Tile3d<double> >(d2, d2, d2);
  std::shared_ptr<Tile3d<double> > forcesy =
    std::make_shared<Tile3d<double> >(d2, d2, d2);
  std::shared_ptr<Tile3d<double> > forcesz =
    std::make_shared<Tile3d<double> >(d2, d2, d2);

  for (i = 0; i < d2; ++i) {
    for (j = 0; j < d2; ++j) {
      for (k = 0; k < d2; ++k) {
        (*forcesx)(i,j,k) = 0.0;
        (*forcesy)(i,j,k) = 0.0;
        (*forcesz)(i,j,k) = 0.0;
      }
    }
  }

  for(i=0; i<block_dim; i++) {
    for(j=0; j<block_dim; j++) {
      for(k=0; k<block_dim; k++) {

        double  x[8],  y[8],  z[8] ;
        double  xd[8],  yd[8],  zd[8] ;
        double  hx[8],  hy[8],  hz[8] ;
        double  sx[8],  sy[8],  sz[8] ;

        x[0] = (*posx)(i,j,k);
        x[1] = (*posx)(i,j,k+1);
        x[2] = (*posx)(i,j+1,k+1);
        x[3] = (*posx)(i,j+1,k);
        x[4] = (*posx)(i+1,j,k);
        x[5] = (*posx)(i+1,j,k+1);
        x[6] = (*posx)(i+1,j+1,k+1);
        x[7] = (*posx)(i+1,j+1,k);

        y[0] = (*posy)(i,j,k);
        y[1] = (*posy)(i,j,k+1);
        y[2] = (*posy)(i,j+1,k+1);
        y[3] = (*posy)(i,j+1,k);
        y[4] = (*posy)(i+1,j,k);
        y[5] = (*posy)(i+1,j,k+1);
        y[6] = (*posy)(i+1,j+1,k+1);
        y[7] = (*posy)(i+1,j+1,k);

        z[0] = (*posz)(i,j,k);
        z[1] = (*posz)(i,j,k+1);
        z[2] = (*posz)(i,j+1,k+1);
        z[3] = (*posz)(i,j+1,k);
        z[4] = (*posz)(i+1,j,k);
        z[5] = (*posz)(i+1,j,k+1);
        z[6] = (*posz)(i+1,j+1,k+1);
        z[7] = (*posz)(i+1,j+1,k);

        xd[0] = (*velx)(i,j,k);
        xd[1] = (*velx)(i,j,k+1);
        xd[2] = (*velx)(i,j+1,k+1);
        xd[3] = (*velx)(i,j+1,k);
        xd[4] = (*velx)(i+1,j,k);
        xd[5] = (*velx)(i+1,j,k+1);
        xd[6] = (*velx)(i+1,j+1,k+1);
        xd[7] = (*velx)(i+1,j+1,k);

        yd[0] = (*vely)(i,j,k);
        yd[1] = (*vely)(i,j,k+1);
        yd[2] = (*vely)(i,j+1,k+1);
        yd[3] = (*vely)(i,j+1,k);
        yd[4] = (*vely)(i+1,j,k);
        yd[5] = (*vely)(i+1,j,k+1);
        yd[6] = (*vely)(i+1,j+1,k+1);
        yd[7] = (*vely)(i+1,j+1,k);

        zd[0] = (*velz)(i,j,k);
        zd[1] = (*velz)(i,j,k+1);
        zd[2] = (*velz)(i,j+1,k+1);
        zd[3] = (*velz)(i,j+1,k);
        zd[4] = (*velz)(i+1,j,k);
        zd[5] = (*velz)(i+1,j,k+1);
        zd[6] = (*velz)(i+1,j+1,k+1);
        zd[7] = (*velz)(i+1,j+1,k);

        double stress = - ( (*pressure_blk)(i,j,k) + (*viscosity_blk)(i,j,k) );

        stress_partial_work(stress, &x[0], &y[0], &z[0], &sx[0], &sy[0], &sz[0]);

        double determ = c.mydomain->element_initial_volume[i][j][k]
                        * (*volume_blk)(i,j,k);
        double precoef = -3.0 * 0.01 * (*ss_blk)(i,j,k) * c.mydomain->element_mass[i][j][k];

        hourglass_partial_work(&x[0], &y[0], &z[0], &xd[0], &yd[0], &zd[0], precoef, determ,
                               &hx[0], &hy[0], &hz[0]);

        (*forcesx)(i,j,k) += hx[0] + sx[0];
        (*forcesy)(i,j,k) += hy[0] + sy[0];
        (*forcesz)(i,j,k) += hz[0] + sz[0];
        (*forcesx)(i,j,k+1) += hx[1] + sx[1];
        (*forcesy)(i,j,k+1) += hy[1] + sy[1];
        (*forcesz)(i,j,k+1) += hz[1] + sz[1];
        (*forcesx)(i,j+1,k+1) += hx[2] + sx[2];
        (*forcesy)(i,j+1,k+1) += hy[2] + sy[2];
        (*forcesz)(i,j+1,k+1) += hz[2] + sz[2];
        (*forcesx)(i,j+1,k) += hx[3] + sx[3];
        (*forcesy)(i,j+1,k) += hy[3] + sy[3];
        (*forcesz)(i,j+1,k) += hz[3] + sz[3];
        (*forcesx)(i+1,j,k) += hx[4] + sx[4];
        (*forcesy)(i+1,j,k) += hy[4] + sy[4];
        (*forcesz)(i+1,j,k) += hz[4] + sz[4];
        (*forcesx)(i+1,j,k+1) += hx[5] + sx[5];
        (*forcesy)(i+1,j,k+1) += hy[5] + sy[5];
        (*forcesz)(i+1,j,k+1) += hz[5] + sz[5];
        (*forcesx)(i+1,j+1,k+1) += hx[6] + sx[6];
        (*forcesy)(i+1,j+1,k+1) += hy[6] + sy[6];
        (*forcesz)(i+1,j+1,k+1) += hz[6] + sz[6];
        (*forcesx)(i+1,j+1,k) += hx[7] + sx[7];
        (*forcesy)(i+1,j+1,k) += hy[7] + sy[7];
        (*forcesz)(i+1,j+1,k) += hz[7] + sz[7];
      }
    }
  }
  c.forcex.put(Quad(q[0], q[1], q[2], q[3]), forcesx);
  c.forcey.put(Quad(q[0], q[1], q[2], q[3]), forcesy);
  c.forcez.put(Quad(q[0], q[1], q[2], q[3]), forcesz);

  return CnC::CNC_Success;
}

template< class dependency_consumer >
void force_reduction_tuner::depends( const Quad & q, lulesh_context & c, dependency_consumer & dC ) const {

  dC.depends(c.forcex, Quad(q[0], q[1], q[2], q[3]));
  dC.depends(c.forcey, Quad(q[0], q[1], q[2], q[3]));
  dC.depends(c.forcez, Quad(q[0], q[1], q[2], q[3]));

  if (q[1]) {
    dC.depends(c.forcex, Quad(q[0], q[1]-1, q[2], q[3]));
    dC.depends(c.forcey, Quad(q[0], q[1]-1, q[2], q[3]));
    dC.depends(c.forcez, Quad(q[0], q[1]-1, q[2], q[3]));
  }
  if (q[2]) {
    dC.depends(c.forcex, Quad(q[0], q[1], q[2]-1, q[3]));
    dC.depends(c.forcey, Quad(q[0], q[1], q[2]-1, q[3]));
    dC.depends(c.forcez, Quad(q[0], q[1], q[2]-1, q[3]));
  }
  if (q[3]) {
    dC.depends(c.forcex, Quad(q[0], q[1], q[2], q[3]-1));
    dC.depends(c.forcey, Quad(q[0], q[1], q[2], q[3]-1));
    dC.depends(c.forcez, Quad(q[0], q[1], q[2], q[3]-1));
  }
}

// Computes node data for: Node ID based tags
int compute_node_data::execute(const Quad &q, lulesh_context & c) const {
  const int iteration = q[0];
  const int block_dim = c.block_dim;
  const int dim = block_dim+1;  // blocksize+1 because nodes
  int xlow, xhigh, ylow, yhigh, zlow, zhigh;
  int i,j,k;
  int nzx, nzy, nzz, cases;
  double delta_time;

  xlow = block_dim*q[1];
  ylow = block_dim*q[2];
  zlow = block_dim*q[3];
  // Add one to handle border cases
  xhigh = xlow + block_dim + (q[1] == c.maxblockid ? 1 : 0);
  yhigh = ylow + block_dim + (q[2] == c.maxblockid ? 1 : 0);
  zhigh = zlow + block_dim + (q[3] == c.maxblockid ? 1 : 0);

  std::shared_ptr<Tile3d<double> > fx, fy, fz;
  std::shared_ptr<Tile3d<double> > posx, posy, posz, velx, vely, velz;

  // Get my force block
  c.forcex.get(Quad(q[0], q[1], q[2], q[3]), fx);
  c.forcey.get(Quad(q[0], q[1], q[2], q[3]), fy);
  c.forcez.get(Quad(q[0], q[1], q[2], q[3]), fz);

  // Store final reduced force here
  double fxfinal[dim][dim][dim];
  double fyfinal[dim][dim][dim];
  double fzfinal[dim][dim][dim];

  for (i = 0; i < dim; ++i) {
    for (j = 0; j < dim; ++j) {
      for (k = 0; k < dim; ++k) {
        fxfinal[i][j][k] = (*fx)(i,j,k);
        fyfinal[i][j][k] = (*fy)(i,j,k);
        fzfinal[i][j][k] = (*fz)(i,j,k);
      }
    }
  }

  // Need force neighbors from -xyz dimensions
  if (q[1]) {
    std::shared_ptr<Tile3d<double> > fxtrans, fytrans, fztrans;
    c.forcex.get(Quad(q[0], q[1]-1, q[2], q[3]), fxtrans);
    c.forcey.get(Quad(q[0], q[1]-1, q[2], q[3]), fytrans);
    c.forcez.get(Quad(q[0], q[1]-1, q[2], q[3]), fztrans);
    for (j = 0; j < dim; ++j) {
      for (k = 0; k < dim; ++k) {
        fxfinal[0][j][k] += (*fxtrans)(dim-1, j, k);
        fyfinal[0][j][k] += (*fytrans)(dim-1, j, k);
        fzfinal[0][j][k] += (*fztrans)(dim-1, j, k);
      }
    }
  }
  if (q[2]) {
    std::shared_ptr<Tile3d<double> > fxtrans, fytrans, fztrans;
    c.forcex.get(Quad(q[0], q[1], q[2]-1, q[3]), fxtrans);
    c.forcey.get(Quad(q[0], q[1], q[2]-1, q[3]), fytrans);
    c.forcez.get(Quad(q[0], q[1], q[2]-1, q[3]), fztrans);
    for (i = 0; i < dim; ++i) {
      for (k = 0; k < dim; ++k) {
        fxfinal[i][0][k] += (*fxtrans)(i, dim-1, k);
        fyfinal[i][0][k] += (*fytrans)(i, dim-1, k);
        fzfinal[i][0][k] += (*fztrans)(i, dim-1, k);
      }
    }
  }
  if (q[3]) {
    std::shared_ptr<Tile3d<double> > fxtrans, fytrans, fztrans;
    c.forcex.get(Quad(q[0], q[1], q[2], q[3]-1), fxtrans);
    c.forcey.get(Quad(q[0], q[1], q[2], q[3]-1), fytrans);
    c.forcez.get(Quad(q[0], q[1], q[2], q[3]-1), fztrans);
    for (i = 0; i < dim; ++i) {
      for (j = 0; j < dim; ++j) {
        fxfinal[i][j][0] += (*fxtrans)(i, j, dim-1);
        fyfinal[i][j][0] += (*fytrans)(i, j, dim-1);
        fzfinal[i][j][0] += (*fztrans)(i, j, dim-1);
      }
    }
  }
  if (q[1] && q[2]) {
    std::shared_ptr<Tile3d<double> > fxtrans, fytrans, fztrans;
    c.forcex.get(Quad(q[0], q[1]-1, q[2]-1, q[3]), fxtrans);
    c.forcey.get(Quad(q[0], q[1]-1, q[2]-1, q[3]), fytrans);
    c.forcez.get(Quad(q[0], q[1]-1, q[2]-1, q[3]), fztrans);
    for (k = 0; k < dim; ++k) {
      fxfinal[0][0][k] += (*fxtrans)(dim-1, dim-1, k);
      fyfinal[0][0][k] += (*fytrans)(dim-1, dim-1, k);
      fzfinal[0][0][k] += (*fztrans)(dim-1, dim-1, k);
    }
  }
  if (q[1] && q[3]) {
    std::shared_ptr<Tile3d<double> > fxtrans, fytrans, fztrans;
    c.forcex.get(Quad(q[0], q[1]-1, q[2], q[3]-1), fxtrans);
    c.forcey.get(Quad(q[0], q[1]-1, q[2], q[3]-1), fytrans);
    c.forcez.get(Quad(q[0], q[1]-1, q[2], q[3]-1), fztrans);
    for (j = 0; j < dim; ++j) {
      fxfinal[0][j][0] += (*fxtrans)(dim-1, j, dim-1);
      fyfinal[0][j][0] += (*fytrans)(dim-1, j, dim-1);
      fzfinal[0][j][0] += (*fztrans)(dim-1, j, dim-1);
    }
  }
  if (q[2] && q[3]) {
    std::shared_ptr<Tile3d<double> > fxtrans, fytrans, fztrans;
    c.forcex.get(Quad(q[0], q[1], q[2]-1, q[3]-1), fxtrans);
    c.forcey.get(Quad(q[0], q[1], q[2]-1, q[3]-1), fytrans);
    c.forcez.get(Quad(q[0], q[1], q[2]-1, q[3]-1), fztrans);
    for (i = 0; i < dim; ++i) {
      fxfinal[i][0][0] += (*fxtrans)(i, dim-1, dim-1);
      fyfinal[i][0][0] += (*fytrans)(i, dim-1, dim-1);
      fzfinal[i][0][0] += (*fztrans)(i, dim-1, dim-1);
    }
  }
  if (q[1] && q[2] && q[3]) {
    std::shared_ptr<Tile3d<double> > fxtrans, fytrans, fztrans;
    c.forcex.get(Quad(q[0], q[1]-1, q[2]-1, q[3]-1), fxtrans);
    c.forcey.get(Quad(q[0], q[1]-1, q[2]-1, q[3]-1), fytrans);
    c.forcez.get(Quad(q[0], q[1]-1, q[2]-1, q[3]-1), fztrans);
    fxfinal[0][0][0] += (*fxtrans)(dim-1, dim-1, dim-1);
    fyfinal[0][0][0] += (*fytrans)(dim-1, dim-1, dim-1);
    fzfinal[0][0][0] += (*fztrans)(dim-1, dim-1, dim-1);
  }

  c.delta_time.get(iteration, delta_time);

  // Get my node tile
  c.posx.get(Quad(q[0]-1, q[1], q[2], q[3]), posx);
  c.posy.get(Quad(q[0]-1, q[1], q[2], q[3]), posy);
  c.posz.get(Quad(q[0]-1, q[1], q[2], q[3]), posz);
  c.velx.get(Quad(q[0]-1, q[1], q[2], q[3]), velx);
  c.vely.get(Quad(q[0]-1, q[1], q[2], q[3]), vely);
  c.velz.get(Quad(q[0]-1, q[1], q[2], q[3]), velz);

  // update acc, velocity, position
  for (i = xlow; i < xhigh; ++i) {
    for (j = ylow; j < yhigh; ++j) {
      for (k = zlow; k < zhigh; ++k) {
        double accx, accy, accz;
        int ii = i - xlow;
        int jj = j - ylow;
        int kk = k - zlow;

        //todo update kinetic values here -- A = F/m (save in F)
        accx = fxfinal[ii][jj][kk] / c.mydomain->node_mass[i][j][k];
        accy = fyfinal[ii][jj][kk] / c.mydomain->node_mass[i][j][k];
        accz = fzfinal[ii][jj][kk] / c.mydomain->node_mass[i][j][k];

        // Acceleration Boundary conditions
        if((i) == 0) accz = 0.0;
        if((j) == 0) accy = 0.0;
        if((k) == 0) accx = 0.0;

        // Velocity calculation
        (*velx)(ii,jj,kk) = (*velx)(ii,jj,kk) + delta_time*accx;
        (*vely)(ii,jj,kk) = (*vely)(ii,jj,kk) + delta_time*accy;
        (*velz)(ii,jj,kk) = (*velz)(ii,jj,kk) + delta_time*accz;

        // Velocity boundary conditions
        if(fabs((*velx)(ii,jj,kk) < 1.0e-7)) (*velx)(ii,jj,kk) = 0.0;
        if(fabs((*vely)(ii,jj,kk) < 1.0e-7)) (*vely)(ii,jj,kk) = 0.0;
        if(fabs((*velz)(ii,jj,kk) < 1.0e-7)) (*velz)(ii,jj,kk) = 0.0;

        // Position Calculation
        (*posx)(ii,jj,kk) += delta_time * (*velx)(ii,jj,kk);
        (*posy)(ii,jj,kk) += delta_time * (*vely)(ii,jj,kk);
        (*posz)(ii,jj,kk) += delta_time * (*velz)(ii,jj,kk);
      }
    }
  }
  c.posx.put(Quad(q[0], q[1], q[2], q[3]), posx);
  c.posy.put(Quad(q[0], q[1], q[2], q[3]), posy);
  c.posz.put(Quad(q[0], q[1], q[2], q[3]), posz);
  c.velx.put(Quad(q[0], q[1], q[2], q[3]), velx);
  c.vely.put(Quad(q[0], q[1], q[2], q[3]), vely);
  c.velz.put(Quad(q[0], q[1], q[2], q[3]), velz);
  return CnC::CNC_Success;
}

template< class dependency_consumer >
void ele1_dep_tuner::depends( const Quad & q, lulesh_context & c, dependency_consumer & dC ) const {
  int d2 = c.block_dim + 1;
  int i,j,k;
  int xlow = c.block_dim*q[1];
  int ylow = c.block_dim*q[2];
  int zlow = c.block_dim*q[3];

  for (i = 0; i < d2; ++i) {
    for (j = 0; j < d2; ++j) {
      for (k = 0; k < d2; ++k) {
        dC.depends(c.posx, Quad(q[0], i+xlow, j+ylow, k+zlow));
        dC.depends(c.posy, Quad(q[0], i+xlow, j+ylow, k+zlow));
        dC.depends(c.posz, Quad(q[0], i+xlow, j+ylow, k+zlow));
        dC.depends(c.velx, Quad(q[0], i+xlow, j+ylow, k+zlow));
        dC.depends(c.vely, Quad(q[0], i+xlow, j+ylow, k+zlow));
        dC.depends(c.velz, Quad(q[0], i+xlow, j+ylow, k+zlow));
      }
    }
  }
}

// Compute element data for: Element ID based tags
int compute_element_step1::execute(const Quad &q, lulesh_context & c) const {
  int iteration = q[0];
  const int block_dim = c.block_dim;
  const int dim = block_dim+1;
  const int maxblockid = c.maxblockid;
  int xlow, xhigh, ylow, yhigh, zlow, zhigh;
  int i,j,k;
  xlow = block_dim*q[1];
  ylow = block_dim*q[2];
  zlow = block_dim*q[3];
  xhigh = xlow + block_dim;
  yhigh = ylow + block_dim;
  zhigh = zlow + block_dim;

  double volume[block_dim*block_dim*block_dim];
  double element_volume_derivative[block_dim*block_dim*block_dim];
  double characteristic_length[block_dim*block_dim*block_dim];
  double delta_time;

  std::shared_ptr<Tile3d<double> > posx, posy, posz, velx, vely, velz;

  // Get my node tile
  c.posx.get(Quad(q[0], q[1], q[2], q[3]), posx);
  c.posy.get(Quad(q[0], q[1], q[2], q[3]), posy);
  c.posz.get(Quad(q[0], q[1], q[2], q[3]), posz);
  c.velx.get(Quad(q[0], q[1], q[2], q[3]), velx);
  c.vely.get(Quad(q[0], q[1], q[2], q[3]), vely);
  c.velz.get(Quad(q[0], q[1], q[2], q[3]), velz);

  // for (i = 0; i < dim; ++i) {
  //   for (j = 0; j < dim; ++j) {
  //     for (k = 0; k < dim; ++k) {
  //       std::cout << "pos: " << (*posx)(i,j,k) << "," << (*posy)(i,j,k) << "," << (*posz)(i,j,k) << "\n";
  //     }
  //   }
  // }
  // update node tile from +xyz neighbors
  if(q[1] < maxblockid) {
    std::shared_ptr<Tile3d<double> > px, py, pz, vx, vy, vz;
    c.posx.get(Quad(q[0], q[1]+1, q[2], q[3]), px);
    c.posy.get(Quad(q[0], q[1]+1, q[2], q[3]), py);
    c.posz.get(Quad(q[0], q[1]+1, q[2], q[3]), pz);
    c.velx.get(Quad(q[0], q[1]+1, q[2], q[3]), vx);
    c.vely.get(Quad(q[0], q[1]+1, q[2], q[3]), vy);
    c.velz.get(Quad(q[0], q[1]+1, q[2], q[3]), vz);

    for (j = 0; j < dim; ++j) {
      for (k = 0; k < dim; ++k) {
        (*posx)(block_dim,j,k) = (*px)(0, j, k);
        (*posy)(block_dim,j,k) = (*py)(0, j, k);
        (*posz)(block_dim,j,k) = (*pz)(0, j, k);
        (*velx)(block_dim,j,k) = (*vx)(0, j, k);
        (*vely)(block_dim,j,k) = (*vy)(0, j, k);
        (*velz)(block_dim,j,k) = (*vz)(0, j, k);
      }
    }
  }
  if(q[2] < maxblockid) {
    std::shared_ptr<Tile3d<double> > px, py, pz, vx, vy, vz;
    c.posx.get(Quad(q[0], q[1], q[2]+1, q[3]), px);
    c.posy.get(Quad(q[0], q[1], q[2]+1, q[3]), py);
    c.posz.get(Quad(q[0], q[1], q[2]+1, q[3]), pz);
    c.velx.get(Quad(q[0], q[1], q[2]+1, q[3]), vx);
    c.vely.get(Quad(q[0], q[1], q[2]+1, q[3]), vy);
    c.velz.get(Quad(q[0], q[1], q[2]+1, q[3]), vz);

    for (i = 0; i < dim; ++i) {
      for (k = 0; k < dim; ++k) {
        (*posx)(i,block_dim,k) = (*px)(i, 0, k);
        (*posy)(i,block_dim,k) = (*py)(i, 0, k);
        (*posz)(i,block_dim,k) = (*pz)(i, 0, k);
        (*velx)(i,block_dim,k) = (*vx)(i, 0, k);
        (*vely)(i,block_dim,k) = (*vy)(i, 0, k);
        (*velz)(i,block_dim,k) = (*vz)(i, 0, k);
      }
    }
  }
  if(q[3] < maxblockid) {
    std::shared_ptr<Tile3d<double> > px, py, pz, vx, vy, vz;
    c.posx.get(Quad(q[0], q[1], q[2], q[3]+1), px);
    c.posy.get(Quad(q[0], q[1], q[2], q[3]+1), py);
    c.posz.get(Quad(q[0], q[1], q[2], q[3]+1), pz);
    c.velx.get(Quad(q[0], q[1], q[2], q[3]+1), vx);
    c.vely.get(Quad(q[0], q[1], q[2], q[3]+1), vy);
    c.velz.get(Quad(q[0], q[1], q[2], q[3]+1), vz);
    for (i = 0; i < dim; ++i) {
      for (j = 0; j < dim; ++j) {
        (*posx)(i,j,block_dim) = (*px)(i, j, 0);
        (*posy)(i,j,block_dim) = (*py)(i, j, 0);
        (*posz)(i,j,block_dim) = (*pz)(i, j, 0);
        (*velx)(i,j,block_dim) = (*vx)(i, j, 0);
        (*vely)(i,j,block_dim) = (*vy)(i, j, 0);
        (*velz)(i,j,block_dim) = (*vz)(i, j, 0);
      }
    }
  }
  if((q[1] < maxblockid) && (q[2] < maxblockid)) {
    std::shared_ptr<Tile3d<double> > px, py, pz, vx, vy, vz;
    c.posx.get(Quad(q[0], q[1]+1, q[2]+1, q[3]), px);
    c.posy.get(Quad(q[0], q[1]+1, q[2]+1, q[3]), py);
    c.posz.get(Quad(q[0], q[1]+1, q[2]+1, q[3]), pz);
    c.velx.get(Quad(q[0], q[1]+1, q[2]+1, q[3]), vx);
    c.vely.get(Quad(q[0], q[1]+1, q[2]+1, q[3]), vy);
    c.velz.get(Quad(q[0], q[1]+1, q[2]+1, q[3]), vz);
    for (k = 0; k < dim; ++k) {
      (*posx)(block_dim,block_dim,k) = (*px)(0, 0, k);
      (*posy)(block_dim,block_dim,k) = (*py)(0, 0, k);
      (*posz)(block_dim,block_dim,k) = (*pz)(0, 0, k);
      (*velx)(block_dim,block_dim,k) = (*vx)(0, 0, k);
      (*vely)(block_dim,block_dim,k) = (*vy)(0, 0, k);
      (*velz)(block_dim,block_dim,k) = (*vz)(0, 0, k);
    }
  }
  if((q[1] < maxblockid) && (q[3] < maxblockid)) {
    std::shared_ptr<Tile3d<double> > px, py, pz, vx, vy, vz;
    c.posx.get(Quad(q[0], q[1]+1, q[2], q[3]+1), px);
    c.posy.get(Quad(q[0], q[1]+1, q[2], q[3]+1), py);
    c.posz.get(Quad(q[0], q[1]+1, q[2], q[3]+1), pz);
    c.velx.get(Quad(q[0], q[1]+1, q[2], q[3]+1), vx);
    c.vely.get(Quad(q[0], q[1]+1, q[2], q[3]+1), vy);
    c.velz.get(Quad(q[0], q[1]+1, q[2], q[3]+1), vz);
    for (j = 0; j < dim; ++j) {
      (*posx)(block_dim,j,block_dim) = (*px)(0, j, 0);
      (*posy)(block_dim,j,block_dim) = (*py)(0, j, 0);
      (*posz)(block_dim,j,block_dim) = (*pz)(0, j, 0);
      (*velx)(block_dim,j,block_dim) = (*vx)(0, j, 0);
      (*vely)(block_dim,j,block_dim) = (*vy)(0, j, 0);
      (*velz)(block_dim,j,block_dim) = (*vz)(0, j, 0);
    }
  }
  if((q[2] < maxblockid) && (q[3] < maxblockid)) {
    std::shared_ptr<Tile3d<double> > px, py, pz, vx, vy, vz;
    c.posx.get(Quad(q[0], q[1], q[2]+1, q[3]+1), px);
    c.posy.get(Quad(q[0], q[1], q[2]+1, q[3]+1), py);
    c.posz.get(Quad(q[0], q[1], q[2]+1, q[3]+1), pz);
    c.velx.get(Quad(q[0], q[1], q[2]+1, q[3]+1), vx);
    c.vely.get(Quad(q[0], q[1], q[2]+1, q[3]+1), vy);
    c.velz.get(Quad(q[0], q[1], q[2]+1, q[3]+1), vz);
    for (i = 0; i < dim; ++i) {
      (*posx)(i,block_dim,block_dim) = (*px)(i, 0, 0);
      (*posy)(i,block_dim,block_dim) = (*py)(i, 0, 0);
      (*posz)(i,block_dim,block_dim) = (*pz)(i, 0, 0);
      (*velx)(i,block_dim,block_dim) = (*vx)(i, 0, 0);
      (*vely)(i,block_dim,block_dim) = (*vy)(i, 0, 0);
      (*velz)(i,block_dim,block_dim) = (*vz)(i, 0, 0);
    }
  }
  if((q[1] < maxblockid) && (q[2] < maxblockid) && (q[3] < maxblockid)) {
    std::shared_ptr<Tile3d<double> > px, py, pz, vx, vy, vz;
    c.posx.get(Quad(q[0], q[1]+1, q[2]+1, q[3]+1), px);
    c.posy.get(Quad(q[0], q[1]+1, q[2]+1, q[3]+1), py);
    c.posz.get(Quad(q[0], q[1]+1, q[2]+1, q[3]+1), pz);
    c.velx.get(Quad(q[0], q[1]+1, q[2]+1, q[3]+1), vx);
    c.vely.get(Quad(q[0], q[1]+1, q[2]+1, q[3]+1), vy);
    c.velz.get(Quad(q[0], q[1]+1, q[2]+1, q[3]+1), vz);
    (*posx)(block_dim,block_dim,block_dim) = (*px)(0, 0, 0);
    (*posy)(block_dim,block_dim,block_dim) = (*py)(0, 0, 0);
    (*posz)(block_dim,block_dim,block_dim) = (*pz)(0, 0, 0);
    (*velx)(block_dim,block_dim,block_dim) = (*vx)(0, 0, 0);
    (*vely)(block_dim,block_dim,block_dim) = (*vy)(0, 0, 0);
    (*velz)(block_dim,block_dim,block_dim) = (*vz)(0, 0, 0);
  }
  c.delta_time.get(iteration, delta_time);

  double p_gradx[block_dim][block_dim][block_dim];
  double p_grady[block_dim][block_dim][block_dim];
  double p_gradz[block_dim][block_dim][block_dim];
  double v_gradx[block_dim][block_dim][block_dim];
  double v_grady[block_dim][block_dim][block_dim];
  double v_gradz[block_dim][block_dim][block_dim];

  std::shared_ptr<Tile3d<double> > volume_block =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > volumeder_block =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > charlen_block =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > pgx =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > pgy =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > pgz =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > vgx =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > vgy =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
  std::shared_ptr<Tile3d<double> > vgz =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);

  double initial_volume = c.mydomain->element_initial_volume[0][0][0];
  for (i = 0; i < block_dim; ++i) {
    for (j = 0; j < block_dim; ++j) {
      for (k = 0; k < block_dim; ++k) {

        double  x[8],  y[8],  z[8] ;
        double  xd[8],  yd[8],  zd[8] ;

        x[0] = (*posx)(i,j,k);
        x[1] = (*posx)(i,j,k+1);
        x[2] = (*posx)(i,j+1,k+1);
        x[3] = (*posx)(i,j+1,k);
        x[4] = (*posx)(i+1,j,k);
        x[5] = (*posx)(i+1,j,k+1);
        x[6] = (*posx)(i+1,j+1,k+1);
        x[7] = (*posx)(i+1,j+1,k);

        y[0] = (*posy)(i,j,k);
        y[1] = (*posy)(i,j,k+1);
        y[2] = (*posy)(i,j+1,k+1);
        y[3] = (*posy)(i,j+1,k);
        y[4] = (*posy)(i+1,j,k);
        y[5] = (*posy)(i+1,j,k+1);
        y[6] = (*posy)(i+1,j+1,k+1);
        y[7] = (*posy)(i+1,j+1,k);

        z[0] = (*posz)(i,j,k);
        z[1] = (*posz)(i,j,k+1);
        z[2] = (*posz)(i,j+1,k+1);
        z[3] = (*posz)(i,j+1,k);
        z[4] = (*posz)(i+1,j,k);
        z[5] = (*posz)(i+1,j,k+1);
        z[6] = (*posz)(i+1,j+1,k+1);
        z[7] = (*posz)(i+1,j+1,k);

        xd[0] = (*velx)(i,j,k);
        xd[1] = (*velx)(i,j,k+1);
        xd[2] = (*velx)(i,j+1,k+1);
        xd[3] = (*velx)(i,j+1,k);
        xd[4] = (*velx)(i+1,j,k);
        xd[5] = (*velx)(i+1,j,k+1);
        xd[6] = (*velx)(i+1,j+1,k+1);
        xd[7] = (*velx)(i+1,j+1,k);

        yd[0] = (*vely)(i,j,k);
        yd[1] = (*vely)(i,j,k+1);
        yd[2] = (*vely)(i,j+1,k+1);
        yd[3] = (*vely)(i,j+1,k);
        yd[4] = (*vely)(i+1,j,k);
        yd[5] = (*vely)(i+1,j,k+1);
        yd[6] = (*vely)(i+1,j+1,k+1);
        yd[7] = (*vely)(i+1,j+1,k);

        zd[0] = (*velz)(i,j,k);
        zd[1] = (*velz)(i,j,k+1);
        zd[2] = (*velz)(i,j+1,k+1);
        zd[3] = (*velz)(i,j+1,k);
        zd[4] = (*velz)(i+1,j,k);
        zd[5] = (*velz)(i+1,j,k+1);
        zd[6] = (*velz)(i+1,j+1,k+1);
        zd[7] = (*velz)(i+1,j+1,k);

        (*volume_block)(i,j,k) = compute_volume_work(x, y, z, initial_volume);
        (*volumeder_block)(i,j,k) = CalcElemDeriv(0.5*delta_time, x, y, z, xd, yd, zd);
        (*charlen_block)(i,j,k) = compute_char_len_work(
                                    x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                                    y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                                    z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7],
                                    (*volume_block)(i,j,k) * initial_volume );
        // std::cout << "volume = " << (*volume_block)(i,j,k) << "\n";
        // std::cout << "volume_derivative = " << (*volumeder_block)(i,j,k) << "\n";
        // std::cout << "characteristic_length = " << (*charlen_block)(i,j,k) << "\n";

        compute_gradients_work(x, y, z, xd, yd, zd, (*volume_block)(i,j,k) * initial_volume,
                               &((*pgx)(i,j,k)), &((*pgy)(i,j,k)), &((*pgz)(i,j,k)), &((*vgx)(i,j,k)), &((*vgy)(i,j,k)), &((*vgz)(i,j,k)));
      }
    }
  }
  c.volume.put(Quad(q[0], q[1], q[2], q[3]), volume_block);
  c.volume_derivative.put(Quad(q[0], q[1], q[2], q[3]), volumeder_block);
  c.characteristic_length.put(Quad(q[0], q[1], q[2], q[3]), charlen_block);
  c.pgx.put(Quad(q[0], q[1], q[2], q[3]), pgx);
  c.pgy.put(Quad(q[0], q[1], q[2], q[3]), pgy);
  c.pgz.put(Quad(q[0], q[1], q[2], q[3]), pgz);
  c.vgx.put(Quad(q[0], q[1], q[2], q[3]), vgx);
  c.vgy.put(Quad(q[0], q[1], q[2], q[3]), vgy);
  c.vgz.put(Quad(q[0], q[1], q[2], q[3]), vgz);

  return CnC::CNC_Success;
}

template< class dependency_consumer >
void ele2_dep_tuner::depends( const Quad & q, lulesh_context & c, dependency_consumer & dC ) const {

  if(q[3] < c.maxblockid)
    dC.depends(c.vgx, Quad(q[0], q[1], q[2], q[3]+1));
  if(q[3] > 0)
    dC.depends(c.vgx, Quad(q[0], q[1], q[2], q[3]-1));
  if(q[2] < c.maxblockid)
    dC.depends(c.vgy, Quad(q[0], q[1], q[2]+1, q[3]));
  if(q[2] > 0)
    dC.depends(c.vgy, Quad(q[0], q[1], q[2]-1, q[3]));
  if(q[1] < c.maxblockid)
    dC.depends(c.vgz, Quad(q[0], q[1]+1, q[2], q[3]));
  if(q[1] > 0)
    dC.depends(c.vgz, Quad(q[0], q[1]-1, q[2], q[3]));
}

// Compute element data for: Element ID based tags
// Only doing element -> element transfers here
int compute_element_step2::execute(const Quad &q, lulesh_context & c) const {
  int iteration = q[0];

  int block_dim = c.block_dim;
  int xlow, xhigh, ylow, yhigh, zlow, zhigh;
  int xvlow, xvhigh, yvlow, yvhigh, zvlow, zvhigh;
  int i,j,k;

  const int dim = block_dim;
  const int d2 = block_dim+2;

  std::shared_ptr<Tile3d<double> >
  energy_block, pressure_block, viscosity_block,
                volume_block, prev_volume_block, volumeder_block,
                charlen_block, pgx, pgy, pgz, vgx, vgy, vgz;

  double vgradx[d2][d2][d2];
  double vgrady[d2][d2][d2];
  double vgradz[d2][d2][d2];

  for (i = 0; i < d2; ++i) {
    for (j = 0; j < d2; ++j) {
      for (k = 0; k < d2; ++k) {
        vgradx[i][j][k] = 0.0;
        vgrady[i][j][k] = 0.0;
        vgradz[i][j][k] = 0.0;
      }
    }
  }

  c.vgx.get(Quad(q[0], q[1], q[2], q[3]), vgx);
  c.vgy.get(Quad(q[0], q[1], q[2], q[3]), vgy);
  c.vgz.get(Quad(q[0], q[1], q[2], q[3]), vgz);
  for (i = 0; i < block_dim; ++i) {
    for (j = 0; j < block_dim; ++j) {
      for (k = 0; k < block_dim; ++k) {
        vgradx[i+1][j+1][k+1] = (*vgx)(i,j,k);
        vgrady[i+1][j+1][k+1] = (*vgy)(i,j,k);
        vgradz[i+1][j+1][k+1] = (*vgz)(i,j,k);
      }
    }
  }

  std::shared_ptr<Tile3d<double> > vxp, vxn, vyp, vyn, vzp, vzn;

  // Get neighbor vgradients
  if(q[3] > 0) {
    c.vgx.get(Quad(q[0], q[1], q[2], q[3]-1), vxn);
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgradx[i+1][j+1][0] = (*vxn)(i,j,dim-1);
      }
    }
  } else {
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgradx[i+1][j+1][0] = vgradx[i+1][j+1][1];
      }
    }
  }
  if(q[3] < c.maxblockid) {
    c.vgx.get(Quad(q[0], q[1], q[2], q[3]+1), vxp);
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgradx[i+1][j+1][d2-1] = (*vxp)(i,j,0);
      }
    }
  }
  if(q[2] > 0) {
    c.vgy.get(Quad(q[0], q[1], q[2]-1, q[3]), vyn);
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgrady[i+1][0][j+1] = (*vyn)(i,dim-1,j);
      }
    }
  } else {
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgrady[i+1][0][j+1] = vgrady[i+1][1][j+1];
      }
    }
  }
  if(q[2] < c.maxblockid) {
    c.vgy.get(Quad(q[0], q[1], q[2]+1, q[3]), vyp);
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgrady[i+1][d2-1][j+1] = (*vyp)(i,0,j);
      }
    }
  }
  if(q[1] > 0) {
    c.vgz.get(Quad(q[0], q[1]-1, q[2], q[3]), vzn);
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgradz[0][i+1][j+1] = (*vzn)(dim-1,i,j);
      }
    }
  } else {
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgradz[0][i+1][j+1] = vgradz[1][i+1][j+1];
      }
    }
  }
  if(q[1] < c.maxblockid) {
    c.vgz.get(Quad(q[0], q[1]+1, q[2], q[3]), vzp);
    for (i = 0; i < block_dim; ++i) {
      for (j = 0; j < block_dim; ++j) {
        vgradz[d2-1][i+1][j+1] = (*vzp)(0,i,j);
      }
    }
  }

  // Get element data
  c.volume.get(Quad(q[0]-1, q[1], q[2], q[3]), prev_volume_block);
  c.volume.get(Quad(q[0], q[1], q[2], q[3]), volume_block);
  c.volume_derivative.get(Quad(q[0], q[1], q[2], q[3]), volumeder_block);
  c.characteristic_length.get(Quad(q[0], q[1], q[2], q[3]), charlen_block);
  c.energy.get(Quad(q[0]-1, q[1], q[2], q[3]), energy_block);
  c.pressure.get(Quad(q[0]-1, q[1], q[2], q[3]), pressure_block);
  c.viscosity.get(Quad(q[0]-1, q[1], q[2], q[3]), viscosity_block);

  c.pgx.get(Quad(q[0], q[1], q[2], q[3]), pgx);
  c.pgy.get(Quad(q[0], q[1], q[2], q[3]), pgy);
  c.pgz.get(Quad(q[0], q[1], q[2], q[3]), pgz);

  double qlin[block_dim][block_dim][block_dim];
  double qquad[block_dim][block_dim][block_dim];

  std::shared_ptr<Tile3d<double> > ss_block =
    std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);

  double minhydro = 1e20;
  double mincour = 1e20;

  for (i = 0; i < block_dim; ++i) {
    for (j = 0; j < block_dim; ++j) {
      for (k = 0; k < block_dim; ++k) {
        double volume, prev_volume, volumeder, charlen, energy, viscosity, pressure, ss;
        double neighborgrads[6];
        double hydro, cour;
        volume = (*volume_block)(i,j,k);
        prev_volume = (*prev_volume_block)(i,j,k);
        volumeder = (*volumeder_block)(i,j,k);
        energy = (*energy_block)(i,j,k);
        viscosity = (*viscosity_block)(i,j,k);
        pressure = (*pressure_block)(i,j,k);
        charlen = (*charlen_block)(i,j,k);

        qlin[i][j][k] = 0.0;
        qquad[i][j][k] = 0.0;

        neighborgrads[0] = vgradx[i+1][j+1][k  ];
        neighborgrads[1] = vgradx[i+1][j+1][k+2];
        neighborgrads[2] = vgrady[i+1][j  ][k+1];
        neighborgrads[3] = vgrady[i+1][j+2][k+1];
        neighborgrads[4] = vgradz[i  ][j+1][k+1];
        neighborgrads[5] = vgradz[i+2][j+1][k+1];

        compute_viscosity_terms_work(volume, volumeder, &neighborgrads[0],
                                     (*pgx)(i,j,k), (*pgy)(i,j,k), (*pgz)(i,j,k),
                                     (*vgx)(i,j,k), (*vgy)(i,j,k), (*vgz)(i,j,k),
                                     &qlin[i][j][k], &qquad[i][j][k]);
        // std::cout << "qquad, qlin = " << qquad[i][j][k] << ", " << qlin[i][j][k] << "\n";

        compute_energy_work(volume, prev_volume, qlin[i][j][k], qquad[i][j][k],
                            &((*energy_block)(i,j,k)), &((*pressure_block)(i,j,k)), &((*viscosity_block)(i,j,k)), &((*ss_block)(i,j,k)));

        compute_time_constraints_work((*ss_block)(i,j,k), volumeder,
                                      charlen, &cour, &hydro);
        mincour = min(mincour, cour);
        minhydro = min(minhydro, hydro);
      }
    }
  }
  if (!q[1] && !q[2] && !q[3]) {
    if (c.debug && (q[0] != c.max_iter))
      std::cout << "Origin energy @ iteration " << iteration << " = " << (*energy_block)(0,0,0) << "\n";
    if (q[0] == c.max_iter)
      std::cout << "Final Origin energy @ iteration " << iteration << " = " << (*energy_block)(0,0,0) << "\n";
  }

  c.sound_speed.put(Quad(q[0], q[1], q[2], q[3]), ss_block);
  c.energy.put(Quad(q[0], q[1], q[2], q[3]), energy_block);
  c.pressure.put(Quad(q[0], q[1], q[2], q[3]), pressure_block);
  c.viscosity.put(Quad(q[0], q[1], q[2], q[3]), viscosity_block);
  c.minhydro.put(Quad(q[0], q[1], q[2], q[3]), minhydro);
  c.mincourant.put(Quad(q[0], q[1], q[2], q[3]), mincour);
  return CnC::CNC_Success;
}

void parse_args(int argc, char **argv, int *p_max_iters, int *p_mesh_size, int *p_block_dim, int *p_num_procs, bool *p_debug) {
  int i;
  *p_max_iters = 10;
  *p_mesh_size = 0;
  *p_block_dim = 0;
  *p_num_procs = 1;
  *p_debug = false;

  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {
      char *argptr;
      *p_max_iters = strtol(argv[++i], &argptr, 10);
    }
    if (strcmp(argv[i], "-n") == 0) {
      char *argptr;
      *p_mesh_size = strtol(argv[++i], &argptr, 10);
    }
    if (strcmp(argv[i], "-b") == 0) {
      char *argptr;
      *p_block_dim = strtol(argv[++i], &argptr, 10);
    }
    if (strcmp(argv[i], "-p") == 0) {
      char *argptr;
      *p_num_procs = strtol(argv[++i], &argptr, 10);
    }
    if (strcmp(argv[i], "-debug") == 0)
      *p_debug = true;
    if (strcmp(argv[i], "-h") == 0)
      printf("usage : %s -n <mesh_size> -b <block_dim> [-p numprocs][-i max_iterations][-debug]\n", argv[0]);
  }
  if( *p_mesh_size == 0 || *p_block_dim == 0 ) {
    printf("usage : %s -n <mesh_size> -b <block_dim> [-p numprocs][-i max_iterations][-debug]\n", argv[0]);
    exit( 11 );
  }
  assert((*p_mesh_size)%(*p_block_dim) == 0);
}

void init_mesh(const int mesh_size, const int block_dim, domain &d, mesh &m) {
  // Initialize the mesh (regular cubes, delta spaced)
  int edge_elements = m.edge_elements = mesh_size;
  int edge_nodes = m.edge_nodes = mesh_size + 1;

  int row_id, column_id, plane_id;
  int i, j, k, element_id, node_id = 0;

  double scale = ((double)edge_elements)/45.0;
  double einit = 3.948746e+7*scale*scale*scale;

  // Initialize the structure; keeps track of size and neighbors
  m.number_nodes = edge_nodes * edge_nodes * edge_nodes;
  m.number_elements = edge_elements * edge_elements * edge_elements;

  m.posx = (double ***) malloc(sizeof(double**)*edge_nodes);
  m.posy = (double ***) malloc(sizeof(double**)*edge_nodes);
  m.posz = (double ***) malloc(sizeof(double**)*edge_nodes);
  d.node_mass = (double ***) malloc(sizeof(double**)*edge_nodes);
  for(i = 0; i < edge_nodes; i++) {
    m.posx[i] = (double **) malloc(sizeof(double*)*edge_nodes);
    m.posy[i] = (double **) malloc(sizeof(double*)*edge_nodes);
    m.posz[i] = (double **) malloc(sizeof(double*)*edge_nodes);
    d.node_mass[i] = (double **) malloc(sizeof(double*)*edge_nodes);
    for(j = 0; j < edge_nodes; j++) {
      m.posx[i][j] = (double *) malloc(sizeof(double)*edge_nodes);
      m.posy[i][j] = (double *) malloc(sizeof(double)*edge_nodes);
      m.posz[i][j] = (double *) malloc(sizeof(double)*edge_nodes);
      d.node_mass[i][j] = (double *) calloc(edge_nodes, sizeof(double));
    }
  }

  d.element_mass = (double ***) malloc(sizeof(double**)*edge_elements);
  d.element_initial_volume = (double ***) malloc(sizeof(double**)*edge_elements);
  for(i = 0; i < edge_elements; i++) {
    d.element_mass[i] = (double **) malloc(sizeof(double*)*edge_elements);
    d.element_initial_volume[i] = (double **) malloc(sizeof(double*)*edge_elements);
    for(j = 0; j < edge_elements; j++) {
      d.element_mass[i][j] = (double *) malloc(sizeof(double)*edge_elements);
      d.element_initial_volume[i][j] = (double *) malloc(sizeof(double)*edge_elements);
    }
  }

  // Setup initial vertices and nodes neighboring nodes
  double delta = 1.125 / ((double) edge_elements);
  for (plane_id = 0; plane_id < edge_nodes; ++plane_id) {
    double z = delta * plane_id;
    for (row_id = 0; row_id < edge_nodes; ++row_id) {
      double y = delta * row_id;
      for (column_id = 0; column_id < edge_nodes; ++column_id) {
        double x = delta * column_id;
        m.posx[plane_id][row_id][column_id] = x;
        m.posy[plane_id][row_id][column_id] = y;
        m.posz[plane_id][row_id][column_id] = z;
      }
    }
  }

  // Setup elements node and element neighbors !!!
  element_id = 0;
  for (plane_id = 0; plane_id < edge_elements; ++plane_id) {
    for (row_id = 0; row_id < edge_elements; ++row_id) {
      for (column_id = 0; column_id < edge_elements; ++column_id) {
        double x[8], y[8], z[8];
        x[0] = m.posx[plane_id][row_id][column_id];
        x[1] = m.posx[plane_id][row_id][column_id+1];
        x[2] = m.posx[plane_id][row_id+1][column_id+1];
        x[3] = m.posx[plane_id][row_id+1][column_id];
        x[4] = m.posx[plane_id+1][row_id][column_id];
        x[5] = m.posx[plane_id+1][row_id][column_id+1];
        x[6] = m.posx[plane_id+1][row_id+1][column_id+1];
        x[7] = m.posx[plane_id+1][row_id+1][column_id];
        y[0] = m.posy[plane_id][row_id][column_id];
        y[1] = m.posy[plane_id][row_id][column_id+1];
        y[2] = m.posy[plane_id][row_id+1][column_id+1];
        y[3] = m.posy[plane_id][row_id+1][column_id];
        y[4] = m.posy[plane_id+1][row_id][column_id];
        y[5] = m.posy[plane_id+1][row_id][column_id+1];
        y[6] = m.posy[plane_id+1][row_id+1][column_id+1];
        y[7] = m.posy[plane_id+1][row_id+1][column_id];
        z[0] = m.posz[plane_id][row_id][column_id];
        z[1] = m.posz[plane_id][row_id][column_id+1];
        z[2] = m.posz[plane_id][row_id+1][column_id+1];
        z[3] = m.posz[plane_id][row_id+1][column_id];
        z[4] = m.posz[plane_id+1][row_id][column_id];
        z[5] = m.posz[plane_id+1][row_id][column_id+1];
        z[6] = m.posz[plane_id+1][row_id+1][column_id+1];
        z[7] = m.posz[plane_id+1][row_id+1][column_id];

        // Compute
        double volume = CalcElemVolume(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                                       y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                                       z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);

        d.element_initial_volume[plane_id][row_id][column_id] = volume;
        d.element_mass[plane_id][row_id][column_id] = volume;

        double node_v = volume / 8.0;
        d.node_mass[plane_id][row_id][column_id] += node_v;
        d.node_mass[plane_id+1][row_id][column_id] += node_v;
        d.node_mass[plane_id][row_id+1][column_id] += node_v;
        d.node_mass[plane_id+1][row_id+1][column_id] += node_v;
        d.node_mass[plane_id][row_id][column_id+1] += node_v;
        d.node_mass[plane_id+1][row_id][column_id+1] += node_v;
        d.node_mass[plane_id][row_id+1][column_id+1] += node_v;
        d.node_mass[plane_id+1][row_id+1][column_id+1] += node_v;
        element_id++;
      }
    }
  }
}

void init_context(lulesh_context &ctx) {
  const int edge_elements = ctx.mymesh->edge_elements;
  const int block_dim = ctx.block_dim;
  const int node_dim = block_dim+1;
  const int num_blocks_per_dim = edge_elements / block_dim;
  const int maxblockid = ctx.maxblockid;

  int edge_nodes = edge_elements + 1;
  double scale = ((double)edge_elements)/45.0;
  double einit = 3.948746e+7*scale*scale*scale;

  int plane_id, row_id, column_id;
  double delta = 1.125 / ((double) edge_elements);

  for (plane_id = 0; plane_id < num_blocks_per_dim; ++plane_id) {
    for (row_id = 0; row_id < num_blocks_per_dim; ++row_id) {
      for (column_id = 0; column_id < num_blocks_per_dim; ++column_id) {
        std::shared_ptr<Tile3d<double> > posx =
          std::make_shared<Tile3d<double> >(node_dim, node_dim, node_dim);
        std::shared_ptr<Tile3d<double> > posy =
          std::make_shared<Tile3d<double> >(node_dim, node_dim, node_dim);
        std::shared_ptr<Tile3d<double> > posz =
          std::make_shared<Tile3d<double> >(node_dim, node_dim, node_dim);
        std::shared_ptr<Tile3d<double> > velx =
          std::make_shared<Tile3d<double> >(node_dim, node_dim, node_dim);
        std::shared_ptr<Tile3d<double> > vely =
          std::make_shared<Tile3d<double> >(node_dim, node_dim, node_dim);
        std::shared_ptr<Tile3d<double> > velz =
          std::make_shared<Tile3d<double> >(node_dim, node_dim, node_dim);
        for (int i = 0; i < node_dim; ++i) {
          for (int j = 0; j < node_dim; ++j) {
            for (int k = 0; k < node_dim; ++k) {
              (*posx)(i,j,k) = ctx.mymesh->posx[i+plane_id*block_dim][j+row_id*block_dim][k+column_id*block_dim];
              (*posy)(i,j,k) = ctx.mymesh->posy[i+plane_id*block_dim][j+row_id*block_dim][k+column_id*block_dim];
              (*posz)(i,j,k) = ctx.mymesh->posz[i+plane_id*block_dim][j+row_id*block_dim][k+column_id*block_dim];
              // std::cout << "pos: " << (*posx)(i,j,k) << "," << (*posy)(i,j,k) << "," << (*posz)(i,j,k) << "\n";
              (*velx)(i,j,k) = 0.0;
              (*vely)(i,j,k) = 0.0;
              (*velz)(i,j,k) = 0.0;
            }
          }
        }
        ctx.posx.put(Quad(0, plane_id, row_id, column_id), posx);
        ctx.posy.put(Quad(0, plane_id, row_id, column_id), posy);
        ctx.posz.put(Quad(0, plane_id, row_id, column_id), posz);
        ctx.velx.put(Quad(0, plane_id, row_id, column_id), velx);
        ctx.vely.put(Quad(0, plane_id, row_id, column_id), vely);
        ctx.velz.put(Quad(0, plane_id, row_id, column_id), velz);
      }
    }
  }

  for (plane_id = 0; plane_id < num_blocks_per_dim; ++plane_id) {
    for (row_id = 0; row_id < num_blocks_per_dim; ++row_id) {
      for (column_id = 0; column_id < num_blocks_per_dim; ++column_id) {
        std::shared_ptr<Tile3d<double> > volume_block =
          std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
        std::shared_ptr<Tile3d<double> > viscosity_block =
          std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
        std::shared_ptr<Tile3d<double> > pressure_block =
          std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
        std::shared_ptr<Tile3d<double> > energy_block =
          std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
        std::shared_ptr<Tile3d<double> > ss_block =
          std::make_shared<Tile3d<double> >(block_dim, block_dim, block_dim);
        for (int i = 0; i < block_dim; ++i) {
          for (int j = 0; j < block_dim; ++j) {
            for (int k = 0; k < block_dim; ++k) {
              (*volume_block)(i,j,k) = 1.0;
              (*viscosity_block)(i,j,k) = 0.0;
              (*pressure_block)(i,j,k) = 0.0;
              (*energy_block)(i,j,k) = (plane_id | row_id | column_id | i | j | k) ? 0.0 : einit;
              (*ss_block)(i,j,k) = 0.0;
            }
          }
        }
        ctx.volume.put(Quad(0, plane_id, row_id, column_id), volume_block);
        ctx.viscosity.put(Quad(0, plane_id, row_id, column_id), viscosity_block);
        ctx.pressure.put(Quad(0, plane_id, row_id, column_id), pressure_block);
        ctx.energy.put(Quad(0, plane_id, row_id, column_id), energy_block);
        ctx.sound_speed.put(Quad(0, plane_id, row_id, column_id), ss_block);
      }
    }
  }

  for(int i = 0; i < edge_nodes; i++) {
    for(int j = 0; j < edge_nodes; j++) {
      free(ctx.mymesh->posx[i][j]);
      free(ctx.mymesh->posy[i][j]);
      free(ctx.mymesh->posz[i][j]);
    }
    free(ctx.mymesh->posx[i]);
    free(ctx.mymesh->posy[i]);
    free(ctx.mymesh->posz[i]);
  }
  free(ctx.mymesh->posx);
  free(ctx.mymesh->posy);
  free(ctx.mymesh->posz);

  double delta_time = 0.5*cbrt(ctx.mydomain->element_initial_volume[0][0][0])/sqrt(2.0*einit);
  double elapsed_time = 0.0;
  double rtime = gettime();

  ctx.elapsed_time.put(0, elapsed_time);
  ctx.delta_time.put(0, delta_time);
  ctx.clock_time.put(0, rtime);
}

int main( int argc, char* argv[] ) {

  int ret;
  int i, j, k;
  double t_start, t_end;

  // Problem configuration
  int    max_iterations = 1;
  int    mesh_size = 0;
  int    block_dim = 0;
  int    use_num_procs = 1;
  bool   debug;

  parse_args(argc, argv, &max_iterations, &mesh_size, &block_dim, &use_num_procs, &debug);

  domain d;
  mesh m;

  init_mesh(mesh_size, block_dim, d, m);

  if(use_num_procs) {
    std::cout << "*** Initialize ::: Setting num_threads = " << use_num_procs << "\n";
    std::cout << "*** Initialize ::: Mesh size = " << mesh_size << "^3\n";
    std::cout << "*** Initialize ::: Run for " << max_iterations << " iterations\n";
    CnC::debug::set_num_threads(use_num_procs); //set number of processors to use
  }

  lulesh_context c(max_iterations, block_dim, debug, &d, &m);
  // CnC::debug::trace(c.posx);

  init_context(c);

  c.iteration.put(0);
  c.wait();



  for(int i = 0; i < c.mymesh->edge_nodes; i++) {
    for(int j = 0; j < c.mymesh->edge_nodes; j++) {
      free(c.mydomain->node_mass[i][j]);
    }
    free(c.mydomain->node_mass[i]);
  }
  free(c.mydomain->node_mass);


  for(i = 0; i < c.mymesh->edge_elements; i++) {
    for(j = 0; j < c.mymesh->edge_elements; j++) {
      free(c.mydomain->element_mass[i][j]);
      free(c.mydomain->element_initial_volume[i][j]);
    }
    free(c.mydomain->element_mass[i]);
    free(c.mydomain->element_initial_volume[i]);
  }
  free(c.mydomain->element_mass);
  free(c.mydomain->element_initial_volume);
}

/******************************************************************************
 * LULESH work routines
 */

void hourglass_partial_work(double *x, double *y, double *z, double *xd, double *yd, double *zd,
                            double precoef, double determ, double *hgfx, double *hgfy, double *hgfz) {

  double pfx[8], pfy[8], pfz[8] ;
  // calc volumederivative for hourglass
  hourglass_derivative(&pfx[0], &pfy[0], &pfz[0], &x[0], &y[0], &z[0]);

  const int gamma[4][8] = {
    {1, 1, -1, -1, -1, -1, 1, 1},
    {1, -1, -1, 1, -1, 1, 1, -1},
    {1, -1, 1, -1, 1, -1, 1, -1},
    {-1, 1, -1, 1, 1, -1, 1, -1}
  };

  double volinv = double(1.0)/determ;
  double coefficient = precoef / cbrt(determ);

  double hourgam[8][4];
  // double xd[8], yd[8], zd[8];

  int i;
  for(int i1=0; i1<4; ++i1) {
    double hourmodx = 0.0;
    double hourmody = 0.0;
    double hourmodz = 0.0;
    for(i = 0; i<8; ++i) {
      hourmodx += x[i] * (double) gamma[i1][i];
      hourmody += y[i] * (double) gamma[i1][i];
      hourmodz += z[i] * (double) gamma[i1][i];
    }
    for(i = 0; i<8; ++i) {
      hourgam[i][i1] = gamma[i1][i] - volinv * ( pfx[i]*hourmodx + pfy[i]*hourmody + pfz[i]*hourmodz );
    }
  }

  double hxx[4];
  for(i=0; i<4; ++i) {
    hxx[i] = hourgam[0][i] * xd[0] + hourgam[1][i] * xd[1] +
             hourgam[2][i] * xd[2] + hourgam[3][i] * xd[3] +
             hourgam[4][i] * xd[4] + hourgam[5][i] * xd[5] +
             hourgam[6][i] * xd[6] + hourgam[7][i] * xd[7];
  }
  for(i = 0; i < 8; i++) {
    hgfx[i] = coefficient *
              (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
               hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
  }
  for(i = 0; i < 4; i++) {
    hxx[i] = hourgam[0][i] * yd[0] + hourgam[1][i] * yd[1] +
             hourgam[2][i] * yd[2] + hourgam[3][i] * yd[3] +
             hourgam[4][i] * yd[4] + hourgam[5][i] * yd[5] +
             hourgam[6][i] * yd[6] + hourgam[7][i] * yd[7];
  }

  for(i = 0; i < 8; i++) {
    hgfy[i] = coefficient *
              (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
               hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
  }
  for(i = 0; i < 4; i++) {
    hxx[i] = hourgam[0][i] * zd[0] + hourgam[1][i] * zd[1] +
             hourgam[2][i] * zd[2] + hourgam[3][i] * zd[3] +
             hourgam[4][i] * zd[4] + hourgam[5][i] * zd[5] +
             hourgam[6][i] * zd[6] + hourgam[7][i] * zd[7];
  }
  for(i = 0; i < 8; i++) {
    hgfz[i] = coefficient *
              (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
               hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
  }
}

void stress_partial_work(double stress, double *x1, double *y1, double *z1,
                         double *strx, double *stry, double *strz) {

  int i;
  // [CalcElemNodeNormals]
  double b0[8], b1[8], b2[8];
  for (i = 0; i < 8; ++i) {
    b0[i] = b1[i] = b2[i] = 0.0;
  }
  if (stress)
    CalcElemNodeNormals(&b0[0], &b1[0], &b2[0], x1, y1, z1);

  for (i = 0; i < 8; ++i) {
    strx[i] = b0[i] * -stress;
    stry[i] = b1[i] * -stress;
    strz[i] = b2[i] * -stress;
  }
}

double compute_volume_work(double *x, double *y, double *z, double initial_volume) {
  double relative_volume;
  relative_volume =
    CalcElemVolume(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                   y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                   z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
  relative_volume = relative_volume / initial_volume;
  return fabs(relative_volume - 1.0) < 1e-10 ? 1.0 : relative_volume;
}

void compute_gradients_work(double *x, double *y, double *z, double *xvel, double *yvel, double *zvel,
                            double volume, double *p_gradx, double *p_grady, double *p_gradz, double *v_gradx, double *v_grady, double *v_gradz) {
  const double ptiny = 1.e-36;

  double dxi, dxj, dxk, dyi, dyj, dyk, dzi, dzj, dzk;
  double axi, axj, axk, ayi, ayj, ayk, azi, azj, azk;

  const double norm = 1.0 / ( volume + ptiny ) ;

  dxj = -0.25*((x[0]+x[1]+x[5]+x[4]) - (x[3]+x[2]+x[6]+x[7])) ;
  dyj = -0.25*((y[0]+y[1]+y[5]+y[4]) - (y[3]+y[2]+y[6]+y[7])) ;
  dzj = -0.25*((z[0]+z[1]+z[5]+z[4]) - (z[3]+z[2]+z[6]+z[7])) ;

  dxi = 0.25*((x[1]+x[2]+x[6]+x[5]) - (x[0]+x[3]+x[7]+x[4])) ;
  dyi = 0.25*((y[1]+y[2]+y[6]+y[5]) - (y[0]+y[3]+y[7]+y[4])) ;
  dzi = 0.25*((z[1]+z[2]+z[6]+z[5]) - (z[0]+z[3]+z[7]+z[4])) ;

  dxk = 0.25*((x[4]+x[5]+x[6]+x[7]) - (x[0]+x[1]+x[2]+x[3])) ;
  dyk = 0.25*((y[4]+y[5]+y[6]+y[7]) - (y[0]+y[1]+y[2]+y[3])) ;
  dzk = 0.25*((z[4]+z[5]+z[6]+z[7]) - (z[0]+z[1]+z[2]+z[3])) ;

  axi = dyi*dzj - dzi*dyj ;
  ayi = dzi*dxj - dxi*dzj ;
  azi = dxi*dyj - dyi*dxj ;
  axj = dyj*dzk - dzj*dyk ;
  ayj = dzj*dxk - dxj*dzk ;
  azj = dxj*dyk - dyj*dxk ;
  axk = dyk*dzi - dzk*dyi ;
  ayk = dzk*dxi - dxk*dzi ;
  azk = dxk*dyi - dyk*dxi ;

  *p_gradz = volume / sqrt(axi*axi + ayi*ayi + azi*azi + ptiny);
  *p_gradx = volume / sqrt(axj*axj + ayj*ayj + azj*azj + ptiny);
  *p_grady = volume / sqrt(axk*axk + ayk*ayk + azk*azk + ptiny);

  axi *= norm;
  ayi *= norm;
  azi *= norm;
  axj *= norm;
  ayj *= norm;
  azj *= norm;
  axk *= norm;
  ayk *= norm;
  azk *= norm;

  dxi = 0.25*((xvel[1]+xvel[2]+xvel[6]+xvel[5]) - (xvel[0]+xvel[3]+xvel[7]+xvel[4])) ;
  dyi = 0.25*((yvel[1]+yvel[2]+yvel[6]+yvel[5]) - (yvel[0]+yvel[3]+yvel[7]+yvel[4])) ;
  dzi = 0.25*((zvel[1]+zvel[2]+zvel[6]+zvel[5]) - (zvel[0]+zvel[3]+zvel[7]+zvel[4])) ;

  dxj = -0.25*((xvel[0]+xvel[1]+xvel[5]+xvel[4]) - (xvel[3]+xvel[2]+xvel[6]+xvel[7])) ;
  dyj = -0.25*((yvel[0]+yvel[1]+yvel[5]+yvel[4]) - (yvel[3]+yvel[2]+yvel[6]+yvel[7])) ;
  dzj = -0.25*((zvel[0]+zvel[1]+zvel[5]+zvel[4]) - (zvel[3]+zvel[2]+zvel[6]+zvel[7])) ;

  dxk = 0.25*((xvel[4]+xvel[5]+xvel[6]+xvel[7]) - (xvel[0]+xvel[1]+xvel[2]+xvel[3])) ;
  dyk = 0.25*((yvel[4]+yvel[5]+yvel[6]+yvel[7]) - (yvel[0]+yvel[1]+yvel[2]+yvel[3])) ;
  dzk = 0.25*((zvel[4]+zvel[5]+zvel[6]+zvel[7]) - (zvel[0]+zvel[1]+zvel[2]+zvel[3])) ;

  *v_gradx = axj*dxi + ayj*dyi + azj*dzi ;
  *v_grady = axk*dxj + ayk*dyj + azk*dzj ;
  *v_gradz = axi*dxk + ayi*dyk + azi*dzk ;
}

void compute_viscosity_terms_work(double volume, double volumeder, double *neighborgrads,
                                  double p_gradx, double p_grady, double p_gradz, double v_gradx, double v_grady, double v_gradz,
                                  double *qlin, double *qquad) {
  int face,index, vindex;
  const double ptiny = 1.e-36;
  const double monoq_limiter_mult = 2.0;
  const double monoq_max_slope = 1.0;
  const double qlc_monoq = 0.5;
  const double qqc_monoq = 2.0/3.0;
  double normx, normy, normz, phix, phiy, phiz ;
  double delvx, delvy, delvz;

  /*  phixi     */
  normx = 1.0 / (v_gradx + ptiny ) ;
  normy = 1.0 / (v_grady + ptiny ) ;
  normz = 1.0 / (v_gradz + ptiny ) ;

  // x,y,z neighbors
  neighborgrads[0] *= normx;
  neighborgrads[1] *= normx;
  neighborgrads[2] *= normy;
  neighborgrads[3] *= normy;
  neighborgrads[4] *= normz;
  neighborgrads[5] *= normz;

  phix = 0.5 * ( neighborgrads[0] + neighborgrads[1] ) ;
  phiy = 0.5 * ( neighborgrads[2] + neighborgrads[3] ) ;
  phiz = 0.5 * ( neighborgrads[4] + neighborgrads[5] ) ;

  for (face = 0; face < 6; face++) {
    neighborgrads[face] *= monoq_limiter_mult;
  }

  if(neighborgrads[0] < phix) phix = neighborgrads[0];
  if(neighborgrads[1] < phix) phix = neighborgrads[1];
  if(phix < 0.0) phix = 0.0;
  if(phix > monoq_limiter_mult) phix = monoq_max_slope;

  if(neighborgrads[2] < phiy ) phiy = neighborgrads[2];
  if(neighborgrads[3] < phiy ) phiy = neighborgrads[3];
  if(phiy < 0.0) phiy = 0.0;
  if(phiy > monoq_max_slope) phiy = monoq_max_slope;

  if(neighborgrads[4] < phiz) phiz = neighborgrads[4];
  if(neighborgrads[5] < phiz) phiz = neighborgrads[5];
  if(phiz < 0.0) phiz = 0.0;
  if(phiz > monoq_max_slope) phiz = monoq_max_slope;

  delvx = v_gradx * p_gradx;
  delvy = v_grady * p_grady;
  delvz = v_gradz * p_gradz;

  if ( delvx > 0.0 ) delvx = 0.0 ;
  if ( delvy > 0.0 ) delvy = 0.0 ;
  if ( delvz > 0.0 ) delvz = 0.0 ;
  double rho = 1.0 / volume;

  if (volumeder > 0.0) {
    *qlin = 0.0;
    *qquad   = 0.0;
  } else {
    *qlin = -qlc_monoq * rho *
            (  delvx * (1.0 - phix) +
               delvy * (1.0 - phiy) +
               delvz * (1.0 - phiz)  ) ;

    *qquad = qqc_monoq * rho *
             (  delvx*delvx * (1.0 - phix*phix) +
                delvy*delvy * (1.0 - phiy*phiy) +
                delvz*delvz * (1.0 - phiz*phiz)  ) ;
  }
}

void compute_energy_work(double volume, double prev_volume, double qlin, double qquad,
                         double *energy_, double *pressure_, double *viscosity_, double *ss_) {

  const double eosvmax = 1.0e+9;
  const double eosvmin = 1.0e-9;
  const double pmin = 0.0;
  const double emin = -1.0e+15;
  const double rho0 = 1.0;
  const double cutoff7 = 1.0e-7;

  double delv, previous_energy, previous_pressure, previous_viscosity;

  delv = volume - prev_volume;
  previous_energy = *energy_;
  previous_pressure = *pressure_;
  previous_viscosity = *viscosity_;

  if (volume < eosvmin) {
    volume = eosvmin;
  }

  if (volume > eosvmax) {
    volume = eosvmax;
  }

  // [EvalEOSForElems]
  double compression = 1.0 / volume - 1.;
  double vchalf = volume - delv * .5;
  double comp_half_step = 1.0 / vchalf - 1.0;
  double c1s = 2.0/3.0;
  double work = 0.0;
  double pressure, viscosity, q_tilde;

  if (volume <= eosvmin) {
    comp_half_step = compression; // impossible due to calling func?
  }
  if (volume >= eosvmax) { // impossible due to calling func?
    previous_pressure  = 0.0;
    compression  = 0.0;
    comp_half_step = 0.0;
  }

  double energy = previous_energy - 0.5 * delv *
                  (previous_pressure + previous_viscosity) + 0.5 * work;

  if (energy  < emin) {
    energy = emin ;
  }

  // [CalcPressureForElems]
  double bvc = c1s * (comp_half_step + 1.0);
  double p_half_step = bvc * energy ;

  if(fabs(p_half_step) < cutoff7)
    p_half_step = 0.0;
  if(volume >= eosvmax) // impossible condition here?
    p_half_step = 0.0;
  if(p_half_step < pmin)
    p_half_step = pmin;

  double vhalf = 1.0 / (1.0 + comp_half_step);
  if (delv > 0.0) {
    viscosity = 0.0; // = qq_old = ql_old
  } else {
    double ssc = ( c1s * energy + vhalf * vhalf * bvc * p_half_step ) / rho0;
    if (ssc <= .1111111e-36) {
      ssc = .3333333e-18;
    } else {
      ssc = sqrt(ssc);
    }
    viscosity = (ssc*qlin + qquad) ;
  }

  energy = energy + 0.5 * delv * (3.0 * (previous_pressure + previous_viscosity)
                                  - 4.0 * (p_half_step + viscosity));

  energy += 0.5 * work;

  if (fabs(energy) < cutoff7) {
    energy = 0.0;
  }
  if(energy  < emin ) {
    energy = emin ;
  }

  // [CalcPressureForElems]
  bvc = c1s * (compression + 1.0);
  pressure = bvc * energy;

  if(fabs(pressure) < cutoff7)
    pressure = 0.0;
  if(volume >= eosvmax ) // impossible condition here?
    pressure = 0.0;
  if(pressure < pmin)
    pressure = pmin;

  if (delv > 0.0) {
    q_tilde = 0.0;
  } else {
    double ssc = ( c1s * energy + volume * volume * bvc * pressure ) / rho0;
    if (ssc <= .1111111e-36) {
      ssc = .3333333e-18;
    } else {
      ssc = sqrt(ssc) ;
    }
    q_tilde = (ssc*qlin + qquad) ;
  }

  energy = energy - (7.0 * (previous_pressure + previous_viscosity) - 8.0 *
                     (p_half_step + viscosity) + (pressure + q_tilde)) * delv / 6.0;

  if (fabs(energy) < cutoff7) {
    energy = 0.0;
  }
  if (energy  < emin) {
    energy = emin ;
  }

  // [CalcPressureForElems]
  bvc = c1s * (compression + 1.0);
  pressure = bvc * energy;

  if(fabs(pressure) < cutoff7)
    pressure = 0.0;

  if(volume >= eosvmax) // impossible condition here?
    pressure = 0.0;
  if(pressure < pmin)
    pressure = pmin;

  if (delv <= 0.0) {
    double ssc = (c1s * energy + volume * volume * bvc * pressure ) / rho0 ;
    if (ssc <= .1111111e-36) {
      ssc = .3333333e-18;
    } else {
      ssc = sqrt(ssc);
    }
    viscosity = (ssc*qlin + qquad) ;
    if (fabs(viscosity) < cutoff7) viscosity = 0.0;
  }

  double sound_speed = (c1s * energy + volume * volume * bvc * pressure) / rho0;
  if (sound_speed <= .1111111e-36) {
    sound_speed = .3333333e-18;
  } else {
    sound_speed = sqrt(sound_speed);
  }

  *energy_ = energy;
  *pressure_ = pressure;
  *viscosity_ = viscosity;
  *ss_ = sound_speed;
}

void compute_time_constraints_work(double sound_speed, double volumeder,
                                   double ch_len, double *dtcour, double *dthydr) {

  const double qqc = 2.0;
  const double qqc2 = 64.0 * qqc * qqc;
  const double dvovmax = 0.1;
  double dtcourant, dthydro, dtf;

  dtcourant = 1.0e+20;
  dthydro = 1.0e+20;
  dtf = sound_speed * sound_speed;

  if (volumeder < 0.0) {
    dtf = dtf + qqc2 * volumeder * volumeder *
          ch_len * ch_len;
  }

  dtf = sqrt(dtf) ;
  dtf = ch_len / dtf ;

  if (volumeder != 0.0) {
    dtcourant = dtf;
  }

  // [CalcHydroConstraintForElems]
  if (volumeder != 0.0) {
    dthydro = dvovmax / (fabs(volumeder)+1.e-20);
  }

  *dthydr = dthydro;
  *dtcour = dtcourant;
}









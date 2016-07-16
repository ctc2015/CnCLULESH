

#ifndef LULESH_H
#define LULESH_H

#include <cnc/cnc.h>
#include <cnc/debug.h>
#include "cnc_common.h"

#define SIZE 20

#define min(x,y) ((x) < (y) ? (x) : (y))

struct domain;
struct mesh;
struct lulesh_context;

typedef struct domain {
  double *** node_mass;
  double *** element_mass;                          // Remains constant
  double *** element_initial_volume;
} domain;

typedef struct mesh {
  int edge_nodes;
  int edge_elements;
  int number_nodes;
  int number_elements;
  double ***posx;
  double ***posy;
  double ***posz;
} mesh;

// Step classes

struct iteration_tag_gen_step {
  int execute( const int & t, lulesh_context & c ) const;
};

struct compute_element_forces {
  int execute(const Quad &q, lulesh_context & c) const;
};

struct compute_node_data {
  int execute(const Quad &q, lulesh_context & c) const;
};

struct compute_element_step1 {
  int execute(const Quad &q, lulesh_context & c) const;
};

struct compute_element_step2 {
  int execute(const Quad &q, lulesh_context & c) const;
};

struct single_use : public CnC::hashmap_tuner {
  int get_count( const int & tag ) const {
    return 1;
  }
  int get_count(const Quad & tag) const {
    return 1;
  }
};

// Step/Item Tuners

struct double_use : public CnC::hashmap_tuner {
  int get_count(const Quad & tag) const {
    return 2;
  }
};

struct triple_use : public CnC::hashmap_tuner {
  int get_count(const Quad & tag) const {
    return 3;
  }
};

struct force_reduction_tuner : public CnC::step_tuner<>, public CnC::hashmap_tuner {
  force_reduction_tuner(lulesh_context &c, int mesh_size, int block_size)
    :ctx(c), mesh_size(mesh_size), block_size(block_size) {}

  template< class dependency_consumer >
  void depends( const Quad & q, lulesh_context & c, dependency_consumer & dC ) const;

  int get_count(const Quad & q) const {
    int count = 1;
    int maxblockid = mesh_size/block_size - 1;
    if(q[1] < maxblockid);
      count *= 2;
    if(q[2] < maxblockid);
      count *= 2;
    if(q[3] < maxblockid);
      count *= 2;
    return count;
  }

private:
  lulesh_context & ctx;
  int mesh_size;
  int block_size;
};

struct ele1_dep_tuner : public CnC::step_tuner<> {
  template< class dependency_consumer >
  void depends( const Quad & q, lulesh_context & c, dependency_consumer & dC ) const;
};

struct ele2_dep_tuner : public CnC::step_tuner<> {
  template< class dependency_consumer >
  void depends( const Quad & q, lulesh_context & c, dependency_consumer & dC ) const;
};

// Context Class
struct lulesh_context : public CnC::context< lulesh_context > {
  int matrix_size, block_size, num_blocks;
  bool debug;

  // Tuners
  force_reduction_tuner force_tuner;

  CnC::step_collection< iteration_tag_gen_step > step_checkdt_gentags;
  CnC::step_collection< compute_element_forces > step_compute_element_forces;
  CnC::step_collection< compute_node_data, force_reduction_tuner > step_compute_node_data;
  CnC::step_collection< compute_element_step1 > step_compute_element_step1;
  CnC::step_collection< compute_element_step2 > step_compute_element_step2;

  // Item Collections
  CnC::item_collection< int, double > delta_time;
  CnC::item_collection< int, double, single_use > elapsed_time;
  CnC::item_collection< int, double, single_use > clock_time;

  // per node items
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, force_reduction_tuner > forcex;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, force_reduction_tuner > forcey;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, force_reduction_tuner > forcez;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > posx;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > posy;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > posz;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > velx;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > vely;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > velz;

  // per element items
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, triple_use > volume;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, double_use > viscosity;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, double_use > pressure;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, single_use > energy;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, single_use > sound_speed;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, single_use > volume_derivative;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, single_use > characteristic_length;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, single_use > pgx;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, single_use > pgy;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> >, single_use > pgz;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > vgx;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > vgy;
  CnC::item_collection< Quad, std::shared_ptr<Tile3d<double> > > vgz;
  CnC::item_collection< Quad, double, single_use > mincourant;
  CnC::item_collection< Quad, double, single_use > minhydro;

  // Tag Collections
  CnC::tag_collection<int> iteration;
  CnC::tag_collection<Quad> node_tags;
  CnC::tag_collection<Quad> element_tags;

  // Configuration
  int max_iter, block_dim, maxblockid;

  //constraints
  double max_delta_time, stop_time;

  //constants
  double hgcoef, ss4o3, qstop, monoq_max_slope, monoq_limiter_mult;
  double qlc_monoq, qqc_monoq, qqc, eosvmax, eosvmin, pmin, emin, dvovmax, refdens;

  //cutoffs
  double e, p, q, v, u;

  domain* mydomain;
  mesh* mymesh;

  // Context constructor
  lulesh_context(int _max_iter=1, int _block_dim=0, bool _debug=0, domain* d=NULL, mesh* m=NULL)
    : CnC::context<lulesh_context>(),

      // Initialize tuner
      force_tuner( *this, m->edge_elements, block_dim ),
      // force_tuner( *this, m->edge_elements, block_dim ),

      // Initialize step collection
      step_checkdt_gentags( *this, "Generate tags for node and element steps"),
      step_compute_element_forces( *this, "compute partial forces for elements onto nodes"),
      step_compute_node_data( *this, force_tuner, "reduce force and update node data"),
      step_compute_element_step1( *this, "1st update element data using updated node values"),
      step_compute_element_step2( *this, "2nd update element data using updated element values"),

      // Initialize each item collection
      delta_time( *this, "delta time" ),
      elapsed_time( *this, "elapsed time" ),
      clock_time( *this, "wall clock time" ),
      forcex( *this, "nodal x force", force_tuner),
      forcey( *this, "nodal y force", force_tuner),
      forcez( *this, "nodal z force", force_tuner),
      posx( *this, "nodal position x" ),
      posy( *this, "nodal position y" ),
      posz( *this, "nodal position z" ),
      velx( *this, "nodal velocity x" ),
      vely( *this, "nodal velocity y" ),
      velz( *this, "nodal velocity z" ),
      pgx( *this, "position gradient" ),
      pgy( *this, "position gradient" ),
      pgz( *this, "position gradient" ),
      vgx( *this, "velocity gradient" ),
      vgy( *this, "velocity gradient" ),
      vgz( *this, "velocity gradient" ),
      volume( *this, "element volume" ),
      volume_derivative( *this, "element volume derivative" ),
      characteristic_length( *this, "element characteristic length" ),
      // quadratic_viscosity_term( *this, "element quadratic viscosity term" ),
      // linear_viscosity_term( *this, "element linear viscosity term" ),
      sound_speed( *this, "element sound speed" ),
      viscosity( *this, "element viscosity" ),
      pressure( *this, "element pressure" ),
      energy( *this, "element energy" ),
      mincourant( *this, "element mincourant" ),
      minhydro( *this, "element minhydro" ),

      // Initialize each tag collection
      iteration( *this ),
      node_tags( *this ),
      element_tags( *this ),

      //initial values
      max_iter(_max_iter),
      block_dim(_block_dim),
      debug(_debug)
  {

    step_checkdt_gentags.controls(iteration);
    step_checkdt_gentags.controls(node_tags);
    step_checkdt_gentags.controls(element_tags);

    iteration.prescribes(step_checkdt_gentags, *this);
    step_checkdt_gentags.consumes(mincourant);
    step_checkdt_gentags.consumes(minhydro);
    step_checkdt_gentags.consumes(delta_time);
    step_checkdt_gentags.produces(delta_time);
    step_checkdt_gentags.consumes(elapsed_time);
    step_checkdt_gentags.produces(elapsed_time);
    step_checkdt_gentags.consumes(clock_time);
    step_checkdt_gentags.produces(clock_time);

    element_tags.prescribes( step_compute_element_forces, *this );
    element_tags.prescribes( step_compute_element_step1, *this );
    element_tags.prescribes( step_compute_element_step2, *this );

    node_tags.prescribes(step_compute_node_data, *this );

    step_compute_element_forces.consumes(pressure);
    step_compute_element_forces.consumes(viscosity);
    step_compute_element_forces.consumes(sound_speed);
    step_compute_element_forces.consumes(volume);
    step_compute_element_forces.consumes(posx);
    step_compute_element_forces.consumes(posy);
    step_compute_element_forces.consumes(posz);
    step_compute_element_forces.consumes(velx);
    step_compute_element_forces.consumes(vely);
    step_compute_element_forces.consumes(velz);
    step_compute_element_forces.produces(forcex);
    step_compute_element_forces.produces(forcey);
    step_compute_element_forces.produces(forcez);

    step_compute_node_data.consumes(forcex);
    step_compute_node_data.consumes(forcey);
    step_compute_node_data.consumes(forcez);
    step_compute_node_data.consumes(posx);
    step_compute_node_data.consumes(posy);
    step_compute_node_data.consumes(posz);
    step_compute_node_data.consumes(velx);
    step_compute_node_data.consumes(vely);
    step_compute_node_data.consumes(velz);
    step_compute_node_data.produces(posx);
    step_compute_node_data.produces(posy);
    step_compute_node_data.produces(posz);
    step_compute_node_data.produces(velx);
    step_compute_node_data.produces(vely);
    step_compute_node_data.produces(velz);

    step_compute_element_step1.consumes(posx);
    step_compute_element_step1.consumes(posy);
    step_compute_element_step1.consumes(posz);
    step_compute_element_step1.consumes(velx);
    step_compute_element_step1.consumes(vely);
    step_compute_element_step1.consumes(velz);
    step_compute_element_step1.consumes(delta_time);
    step_compute_element_step1.produces(volume);
    step_compute_element_step1.produces(volume_derivative);
    step_compute_element_step1.produces(characteristic_length);
    step_compute_element_step1.produces(pgx);
    step_compute_element_step1.produces(pgy);
    step_compute_element_step1.produces(pgz);
    step_compute_element_step1.produces(vgx);
    step_compute_element_step1.produces(vgy);
    step_compute_element_step1.produces(vgz);

    step_compute_element_step2.consumes(energy);
    step_compute_element_step2.consumes(pressure);
    step_compute_element_step2.consumes(viscosity);
    step_compute_element_step2.consumes(volume);
    step_compute_element_step2.consumes(volume_derivative);
    step_compute_element_step2.consumes(characteristic_length);
    step_compute_element_step2.consumes(pgx);
    step_compute_element_step2.consumes(pgy);
    step_compute_element_step2.consumes(pgz);
    step_compute_element_step2.consumes(vgx);
    step_compute_element_step2.consumes(vgy);
    step_compute_element_step2.consumes(vgz);
    step_compute_element_step2.produces(mincourant);
    step_compute_element_step2.produces(minhydro);
    step_compute_element_step2.produces(energy);
    step_compute_element_step2.produces(pressure);
    step_compute_element_step2.produces(energy);
    step_compute_element_step2.produces(viscosity);
    step_compute_element_step2.produces(sound_speed);

    mydomain = d;
    mymesh = m;
    maxblockid = m->edge_elements / block_dim - 1;
    num_blocks = m->edge_elements / block_dim;

    //constraints
    max_delta_time = 1.0e-2;
    stop_time = 1.0e-2;
    hgcoef = 3.0;
    ss4o3 = 4.0/3.0;
    qstop = 1.0e+12;
    monoq_max_slope = 1.0;
    monoq_limiter_mult = 2.0;
    qlc_monoq = 0.5;
    qqc_monoq = 2.0/3.0;
    qqc = 2.0;
    eosvmax = 1.0e+9;
    eosvmin = 1.0e-9;
    pmin = 0.0;
    emin = -1.0e+15;
    dvovmax = 0.1;
    refdens = 1.0;

    //cutoffs
    e = 1.0e-7;
    p = 1.0e-7;
    q = 1.0e-7;
    v = 1.0e-10;
    u = 1.0e-7;

#if 0 // Debug
    CnC::debug::collect_scheduler_statistics(*this);
    CnC::debug::trace_all( *this );
#endif
  }
  ~lulesh_context() {
  }

};

/******************************************************************************
 * helper functions
 */

void hourglass_partial_work(double *x, double *y, double *z, double *xd, double *yd, double *zd,
                            double precoef, double determ, double *hgfx, double *hgfy, double *hgfz);

void stress_partial_work(double stress, double *x1, double *y1, double *z1,
                         double *strx, double *stry, double *strz);

double compute_volume_work(double *x, double *y, double *z, double initial_volume);

void compute_viscosity_terms_work(double volume, double volumeder, double *neighborgrads,
                                  double p_gradx, double p_grady, double p_gradz, double v_gradx, double v_grady, double v_gradz,
                                  double *qlin, double *qquad);

void compute_energy_work(double volume, double prev_volume, double qlin, double qquad,
                         double *energy_, double *pressure_, double *viscosity_, double *ss_);

void compute_time_constraints_work(double sound_speed, double volumeder,
                                   double ch_len, double *dtcour, double *dthydr);

void compute_gradients_work(double *x, double *y, double *z, double *xvel, double *yvel, double *zvel,
                            double volume, double *p_gradx, double *p_grady, double *p_gradz, double *v_gradx, double *v_grady, double *v_gradz);

static inline
void VoluDer(const double x0, const double x1, const double x2,
             const double x3, const double x4, const double x5,
             const double y0, const double y1, const double y2,
             const double y3, const double y4, const double y5,
             const double z0, const double z1, const double z2,
             const double z3, const double z4, const double z5,
             double* dvdx, double* dvdy, double* dvdz) {
  const double twelfth = 1.0/12.0 ;
  *dvdx =
    (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
    (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
    (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
  *dvdy =
    - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
    (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
    (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);
  *dvdz =
    - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
    (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
    (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);
  *dvdx *= twelfth;
  *dvdy *= twelfth;
  *dvdz *= twelfth;
}

/******************************************/

static inline
double CalcElemDeriv( double const dt2,
                      double const x[],
                      double const y[],
                      double const z[],
                      double const xvel[],
                      double const yvel[],
                      double const zvel[]) {
  const double x0 = x[0] - dt2*xvel[0] ;
  const double x1 = x[1] - dt2*xvel[1] ;
  const double x2 = x[2] - dt2*xvel[2] ;
  const double x3 = x[3] - dt2*xvel[3] ;
  const double x4 = x[4] - dt2*xvel[4] ;
  const double x5 = x[5] - dt2*xvel[5] ;
  const double x6 = x[6] - dt2*xvel[6] ;
  const double x7 = x[7] - dt2*xvel[7] ;

  const double y0 = y[0] - dt2*yvel[0] ;
  const double y1 = y[1] - dt2*yvel[1] ;
  const double y2 = y[2] - dt2*yvel[2] ;
  const double y3 = y[3] - dt2*yvel[3] ;
  const double y4 = y[4] - dt2*yvel[4] ;
  const double y5 = y[5] - dt2*yvel[5] ;
  const double y6 = y[6] - dt2*yvel[6] ;
  const double y7 = y[7] - dt2*yvel[7] ;

  const double z0 = z[0] - dt2*zvel[0] ;
  const double z1 = z[1] - dt2*zvel[1] ;
  const double z2 = z[2] - dt2*zvel[2] ;
  const double z3 = z[3] - dt2*zvel[3] ;
  const double z4 = z[4] - dt2*zvel[4] ;
  const double z5 = z[5] - dt2*zvel[5] ;
  const double z6 = z[6] - dt2*zvel[6] ;
  const double z7 = z[7] - dt2*zvel[7] ;

  double b[3][8];

  double fjxxi, fjxet, fjxze;
  double fjyxi, fjyet, fjyze;
  double fjzxi, fjzet, fjzze;
  double cjxxi, cjxet, cjxze;
  double cjyxi, cjyet, cjyze;
  double cjzxi, cjzet, cjzze;

  fjxxi = 0.125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
  fjxet = 0.125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
  fjxze = 0.125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

  fjyxi = 0.125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
  fjyet = 0.125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
  fjyze = 0.125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

  fjzxi = 0.125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
  fjzet = 0.125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
  fjzze = 0.125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

  /* compute cofactors */
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);
  b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  b[0][1] =      cjxxi  -  cjxet  -  cjxze;
  b[0][2] =      cjxxi  +  cjxet  -  cjxze;
  b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;

  b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  b[1][1] =      cjyxi  -  cjyet  -  cjyze;
  b[1][2] =      cjyxi  +  cjyet  -  cjyze;
  b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;

  b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  b[2][1] =      cjzxi  -  cjzet  -  cjzze;
  b[2][2] =      cjzxi  +  cjzet  -  cjzze;
  b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;

  /* calculate jacobian determinant (volume) */
  const double inv_detJ = 1.0 / (8.0 * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet));
  // std::cout << "inv_detJ = " << inv_detJ << "\n";
  double dxx, dyy, dzz;

  dxx = inv_detJ * ( b[0][0] * (xvel[0]-xvel[6])
                     + b[0][1] * (xvel[1]-xvel[7])
                     + b[0][2] * (xvel[2]-xvel[4])
                     + b[0][3] * (xvel[3]-xvel[5]) );

  dyy = inv_detJ * ( b[1][0] * (yvel[0]-yvel[6])
                     + b[1][1] * (yvel[1]-yvel[7])
                     + b[1][2] * (yvel[2]-yvel[4])
                     + b[1][3] * (yvel[3]-yvel[5]) );

  dzz = inv_detJ * ( b[2][0] * (zvel[0]-zvel[6])
                     + b[2][1] * (zvel[1]-zvel[7])
                     + b[2][2] * (zvel[2]-zvel[4])
                     + b[2][3] * (zvel[3]-zvel[5]) );

  return dxx + dyy + dzz;
}

/******************************************/

static inline
void hourglass_derivative(double *dvdx, double *dvdy, double *dvdz, double *x1, double *y1, double *z1) {

  VoluDer(x1[1], x1[2], x1[3], x1[4], x1[5], x1[7],
          y1[1], y1[2], y1[3], y1[4], y1[5], y1[7],
          z1[1], z1[2], z1[3], z1[4], z1[5], z1[7],
          &dvdx[0], &dvdy[0], &dvdz[0]);
  VoluDer(x1[0], x1[1], x1[2], x1[7], x1[4], x1[6],
          y1[0], y1[1], y1[2], y1[7], y1[4], y1[6],
          z1[0], z1[1], z1[2], z1[7], z1[4], z1[6],
          &dvdx[3], &dvdy[3], &dvdz[3]);
  VoluDer(x1[3], x1[0], x1[1], x1[6], x1[7], x1[5],
          y1[3], y1[0], y1[1], y1[6], y1[7], y1[5],
          z1[3], z1[0], z1[1], z1[6], z1[7], z1[5],
          &dvdx[2], &dvdy[2], &dvdz[2]);
  VoluDer(x1[2], x1[3], x1[0], x1[5], x1[6], x1[4],
          y1[2], y1[3], y1[0], y1[5], y1[6], y1[4],
          z1[2], z1[3], z1[0], z1[5], z1[6], z1[4],
          &dvdx[1], &dvdy[1], &dvdz[1]);
  VoluDer(x1[7], x1[6], x1[5], x1[0], x1[3], x1[1],
          y1[7], y1[6], y1[5], y1[0], y1[3], y1[1],
          z1[7], z1[6], z1[5], z1[0], z1[3], z1[1],
          &dvdx[4], &dvdy[4], &dvdz[4]);
  VoluDer(x1[4], x1[7], x1[6], x1[1], x1[0], x1[2],
          y1[4], y1[7], y1[6], y1[1], y1[0], y1[2],
          z1[4], z1[7], z1[6], z1[1], z1[0], z1[2],
          &dvdx[5], &dvdy[5], &dvdz[5]);
  VoluDer(x1[5], x1[4], x1[7], x1[2], x1[1], x1[3],
          y1[5], y1[4], y1[7], y1[2], y1[1], y1[3],
          z1[5], z1[4], z1[7], z1[2], z1[1], z1[3],
          &dvdx[6], &dvdy[6], &dvdz[6]);
  VoluDer(x1[6], x1[5], x1[4], x1[3], x1[2], x1[0],
          y1[6], y1[5], y1[4], y1[3], y1[2], y1[0],
          z1[6], z1[5], z1[4], z1[3], z1[2], z1[0],
          &dvdx[7], &dvdy[7], &dvdz[7]);
}

/******************************************/

static inline
void SumElemFaceNormal(double *normalX0, double *normalY0, double *normalZ0,
                       double *normalX1, double *normalY1, double *normalZ1,
                       double *normalX2, double *normalY2, double *normalZ2,
                       double *normalX3, double *normalY3, double *normalZ3,
                       const double x0, const double y0, const double z0,
                       const double x1, const double y1, const double z1,
                       const double x2, const double y2, const double z2,
                       const double x3, const double y3, const double z3) {
  double bisectX0 = 0.5 * (x3 + x2 - x1 - x0);
  double bisectY0 = 0.5 * (y3 + y2 - y1 - y0);
  double bisectZ0 = 0.5 * (z3 + z2 - z1 - z0);
  double bisectX1 = 0.5 * (x2 + x1 - x3 - x0);
  double bisectY1 = 0.5 * (y2 + y1 - y3 - y0);
  double bisectZ1 = 0.5 * (z2 + z1 - z3 - z0);
  double areaX = 0.25 * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
  double areaY = 0.25 * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
  double areaZ = 0.25 * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

  *normalX0 += areaX;
  *normalX1 += areaX;
  *normalX2 += areaX;
  *normalX3 += areaX;

  *normalY0 += areaY;
  *normalY1 += areaY;
  *normalY2 += areaY;
  *normalY3 += areaY;

  *normalZ0 += areaZ;
  *normalZ1 += areaZ;
  *normalZ2 += areaZ;
  *normalZ3 += areaZ;
}

/******************************************/

static inline
void CalcElemNodeNormals(double *pfx, double *pfy, double *pfz,
                         const double x[8], const double y[8], const double z[8]) {
  for (int i = 0 ; i < 8 ; ++i) {
    pfx[i] = 0.0;
    pfy[i] = 0.0;
    pfz[i] = 0.0;
  }
  /* evaluate face one: nodes 0, 1, 2, 3 */
  SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                    &pfx[1], &pfy[1], &pfz[1],
                    &pfx[2], &pfy[2], &pfz[2],
                    &pfx[3], &pfy[3], &pfz[3],
                    x[0], y[0], z[0], x[1], y[1], z[1],
                    x[2], y[2], z[2], x[3], y[3], z[3]);
  /* evaluate face two: nodes 0, 4, 5, 1 */
  SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                    &pfx[4], &pfy[4], &pfz[4],
                    &pfx[5], &pfy[5], &pfz[5],
                    &pfx[1], &pfy[1], &pfz[1],
                    x[0], y[0], z[0], x[4], y[4], z[4],
                    x[5], y[5], z[5], x[1], y[1], z[1]);
  /* evaluate face three: nodes 1, 5, 6, 2 */
  SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
                    &pfx[5], &pfy[5], &pfz[5],
                    &pfx[6], &pfy[6], &pfz[6],
                    &pfx[2], &pfy[2], &pfz[2],
                    x[1], y[1], z[1], x[5], y[5], z[5],
                    x[6], y[6], z[6], x[2], y[2], z[2]);
  /* evaluate face four: nodes 2, 6, 7, 3 */
  SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
                    &pfx[6], &pfy[6], &pfz[6],
                    &pfx[7], &pfy[7], &pfz[7],
                    &pfx[3], &pfy[3], &pfz[3],
                    x[2], y[2], z[2], x[6], y[6], z[6],
                    x[7], y[7], z[7], x[3], y[3], z[3]);
  /* evaluate face five: nodes 3, 7, 4, 0 */
  SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
                    &pfx[7], &pfy[7], &pfz[7],
                    &pfx[4], &pfy[4], &pfz[4],
                    &pfx[0], &pfy[0], &pfz[0],
                    x[3], y[3], z[3], x[7], y[7], z[7],
                    x[4], y[4], z[4], x[0], y[0], z[0]);
  /* evaluate face six: nodes 4, 7, 6, 5 */
  SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
                    &pfx[7], &pfy[7], &pfz[7],
                    &pfx[6], &pfy[6], &pfz[6],
                    &pfx[5], &pfy[5], &pfz[5],
                    x[4], y[4], z[4], x[7], y[7], z[7],
                    x[6], y[6], z[6], x[5], y[5], z[5]);
}

/******************************************/

static inline
double CalcElemVolume(
  const double x0, const double x1, const double x2, const double x3,
  const double x4, const double x5, const double x6, const double x7,
  const double y0, const double y1, const double y2, const double y3,
  const double y4, const double y5, const double y6, const double y7,
  const double z0, const double z1, const double z2, const double z3,
  const double z4, const double z5, const double z6, const double z7 ) {
  double twelveth = 1.0/12.0;
  double dx61 = x6 - x1;
  double dy61 = y6 - y1;
  double dz61 = z6 - z1;
  double dx70 = x7 - x0;
  double dy70 = y7 - y0;
  double dz70 = z7 - z0;
  double dx63 = x6 - x3;
  double dy63 = y6 - y3;
  double dz63 = z6 - z3;
  double dx20 = x2 - x0;
  double dy20 = y2 - y0;
  double dz20 = z2 - z0;
  double dx50 = x5 - x0;
  double dy50 = y5 - y0;
  double dz50 = z5 - z0;
  double dx64 = x6 - x4;
  double dy64 = y6 - y4;
  double dz64 = z6 - z4;
  double dx31 = x3 - x1;
  double dy31 = y3 - y1;
  double dz31 = z3 - z1;
  double dx72 = x7 - x2;
  double dy72 = y7 - y2;
  double dz72 = z7 - z2;
  double dx43 = x4 - x3;
  double dy43 = y4 - y3;
  double dz43 = z4 - z3;
  double dx57 = x5 - x7;
  double dy57 = y5 - y7;
  double dz57 = z5 - z7;
  double dx14 = x1 - x4;
  double dy14 = y1 - y4;
  double dz14 = z1 - z4;
  double dx25 = x2 - x5;
  double dy25 = y2 - y5;
  double dz25 = z2 - z5;
#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
   ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

  double volume =
    TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20, dy31 + dy72, dy63, dy20,
                   dz31 + dz72, dz63, dz20) +
    TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70, dy43 + dy57, dy64, dy70,
                   dz43 + dz57, dz64, dz70) +
    TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50, dy14 + dy25, dy61, dy50,
                   dz14 + dz25, dz61, dz50);

#undef TRIPLE_PRODUCT

  volume *= twelveth;

  return volume ;
}

/******************************************/

static inline
double AreaFace( const double x0, const double x1,
                 const double x2, const double x3,
                 const double y0, const double y1,
                 const double y2, const double y3,
                 const double z0, const double z1,
                 const double z2, const double z3) {
  double fx = (x2 - x0) - (x3 - x1);
  double fy = (y2 - y0) - (y3 - y1);
  double fz = (z2 - z0) - (z3 - z1);
  double gx = (x2 - x0) + (x3 - x1);
  double gy = (y2 - y0) + (y3 - y1);
  double gz = (z2 - z0) + (z3 - z1);
  double area =
    (fx * fx + fy * fy + fz * fz) *
    (gx * gx + gy * gy + gz * gz) -
    (fx * gx + fy * gy + fz * gz) *
    (fx * gx + fy * gy + fz * gz);
  return area ;
}

/******************************************/

static inline
double compute_char_len_work(
  const double x0, const double x1, const double x2, const double x3,
  const double x4, const double x5, const double x6, const double x7,
  const double y0, const double y1, const double y2, const double y3,
  const double y4, const double y5, const double y6, const double y7,
  const double z0, const double z1, const double z2, const double z3,
  const double z4, const double z5, const double z6, const double z7,
  const double volume) {

  double length = 0.0;
  double a;
  a = AreaFace(x0,x1,x2,x3,
               y0,y1,y2,y3,
               z0,z1,z2,z3) ;
  length = std::max(a,length) ;

  a = AreaFace(x4,x5,x6,x7,
               y4,y5,y6,y7,
               z4,z5,z6,z7) ;
  length = std::max(a,length) ;

  a = AreaFace(x0,x1,x5,x4,
               y0,y1,y5,y4,
               z0,z1,z5,z4) ;
  length = std::max(a,length) ;

  a = AreaFace(x1,x2,x6,x5,
               y1,y2,y6,y5,
               z1,z2,z6,z5) ;
  length = std::max(a,length) ;

  a = AreaFace(x2,x3,x7,x6,
               y2,y3,y7,y6,
               z2,z3,z7,z6) ;
  length = std::max(a,length) ;

  a = AreaFace(x3,x0,x4,x7,
               y3,y0,y4,y7,
               z3,z0,z4,z7) ;
  length = std::max(a,length) ;

  return (4.0 * volume / sqrt(length));
}

/******************************************/


#endif // LULESH_H

/*
 * Montee capillaire dans un tube axisymetrique, maillage AMR.
 *
 * Compilation :
 * qcc -autolink -disable-dimensions -O2 -Wall jurin_axi_amr.c -o jurin_axi_amr -lm
 *
 * Execution :
 * ./jurin_axi_amr
 * Sorties : hauteur_axi_amr.dat, cellules.dat, vitesse_t015.dat (unites SI).
 */

#include "grid/quadtree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"

// Proprietes inspirees du cas Omega = 1 de Grunding et al. (2020).
// MU_L adaptee au tube cylindrique pour conserver Omega = 1.
// Difference avec Bastien : fluides, gravite et glissement differents.
#define R           5.e-3                       // Rayon du tube (m)
#define HDOM        (16.*R)                     // Hauteur du domaine (m)
#define H0          (4.*R)                      // Hauteur moyenne initiale (m)
#define THETA       (30.*pi/180.)               // Angle de contact (rad)
#define RHO_L       83.1                        // Masse volumique du liquide (kg/m3)
#define RHO_G       (RHO_L/1000.)               // Masse volumique du gaz (kg/m3)
#define MU_L        0.0026519315195571015      // Viscosite dynamique du liquide (Pa.s)
#define MU_G        (MU_L/1000.)                // Viscosite dynamique du gaz (Pa.s)
#define SIGMA       0.04                        // Tension superficielle (N/m)
#define G           4.17                        // Gravite (m/s2)
#define SLIP        (R/5.)                      // Longueur de glissement (m)
#define T_END       0.8                         // Duree simulee (s)
#define DT_OUTPUT   5.e-4                       // Intervalle entre mesures (s)

#ifndef MINLEVEL
#define MINLEVEL    5                           // Niveau de base ; paroi alignee des le niveau 4
#endif
#ifndef MAXLEVEL
#define MAXLEVEL    9                           // 2^(MAXLEVEL-4) cellules par rayon
#endif
#ifndef F_TOL
#define F_TOL       1.e-3                       // Tolerance sur la fraction liquide
#endif
#ifndef U_TOL
#define U_TOL       1.e-3                       // Tolerance sur chaque composante de vitesse (m/s)
#endif
#define T_COMPARE   0.15                        // Instant de comparaison des vitesses (s)
#define SAMPLE_NX   1024                        // Points fixes suivant x (direction axiale)
#define SAMPLE_NY   64                          // Points fixes suivant y (direction radiale)

#define PIN         (RHO_G*G*HDOM)              // Pression imposee a l'entree (Pa)

vector h[];
h.t[top] = contact_angle (THETA);

// y=0 : axe de symetrie ; y=R : paroi mouillante masquee.
u.n[bottom] = dirichlet (0.);
u.t[bottom] = neumann (0.);
u.n[top] = dirichlet (0.);
u.t[top] = navier (0., SLIP);
// Pression sans gradient radial au mur masque ; evite une forme 0/0.
p[top] = neumann (0.);
pf[top] = neumann (0.);

// x=0 : reservoir liquide ; x=HDOM : sortie gaz.
u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (PIN);
pf[left] = dirichlet (PIN);
f[left] = dirichlet (1.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet (0.);
f[right] = dirichlet (0.);

int main ()
{
  size (HDOM);
  origin (0., 0.);
  init_grid (1 << MINLEVEL);

  rho1 = RHO_L;
  rho2 = RHO_G;
  mu1 = MU_L;
  mu2 = MU_G;
  f.sigma = SIGMA;
  f.height = h;
  TOLERANCE = 1.e-6;

  const double hj = 2.*SIGMA*cos(THETA)/((RHO_L - RHO_G)*G*R);
  const double omega = sqrt(128.*SIGMA*cos(THETA)*sq(MU_L)/
                            (pow(RHO_L,3.)*sq(G)*pow(R,5.)));
  const double dxmin = HDOM/(1 << MAXLEVEL);
  fprintf (stderr,
           "AXI AMR | niveaux=%d..%d | dx_min=%.5g m | cellules/R=%g | "
           "slip/dx_min=%.3g | hJ=%.7g m | Omega=%.6g | F_TOL=%g | U_TOL=%g m/s\n",
           MINLEVEL, MAXLEVEL, dxmin, R/dxmin, SLIP/dxmin, hj, omega,
           F_TOL, U_TOL);
  FILE * fp = fopen ("parametres.json", "w");
  if (!fp) { perror ("parametres.json"); return 1; }
  fprintf (fp,
           "{\"minlevel\":%d,\"maxlevel\":%d,\"f_tol\":%.17g,\"u_tol\":%.17g,"
           "\"t_end\":%.17g,\"t_compare\":%.17g,\"length_m\":%.17g,"
           "\"radius_m\":%.17g,\"sample_nx\":%d,\"sample_ny\":%d}\n",
           MINLEVEL, MAXLEVEL, F_TOL, U_TOL, T_END, T_COMPARE, HDOM, R,
           SAMPLE_NX, SAMPLE_NY);
  fclose (fp);
  run();
}

// Difference avec Bastien : hauteur initiale de 20 mm et menisque courbe.
event init (i = 0)
{
  // Difference avec Bastien : quadtree et paroi definie par mask().
  // Le masque est applique apres la recreation de la grille par run().
  mask (y > R ? top : none);

  const double rc = R/cos(THETA);
  const double q = sqrt(sq(rc) - sq(R));
  const double shift = rc - 2.*(pow(rc,3.) - pow(q,3.))/(3.*sq(R));
  // Prolongement au mur pour les sommets evalues hors de la zone masquee.
  fraction (f, H0 - shift + rc - sqrt(sq(rc) - sq(min(y, R))) - x);

  foreach() {
    u.x[] = 0.;
    u.y[] = 0.;
    p[] = f[]*(PIN - RHO_L*G*x) +
          (1. - f[])*RHO_G*G*(HDOM - x);
  }
  boundary ({f, u, p});
}

event acceleration (i++)
{
  face vector av = a;
  foreach_face(x)
    av.x[] -= G;
}

// Raffinement selon l'interface et la vitesse ; mise a jour des hauteurs de contact.
event adapt (i++)
{
  adapt_wavelet ((scalar *){f, u.x, u.y},
                 (double[]){F_TOL, U_TOL, U_TOL}, MAXLEVEL, MINLEVEL);
  heights (f, h);
}

// Hauteur moyenne par le volume ; hauteur sur l'axe estimee au premier rang.
// Bastien repere plutot la derniere cellule liquide sur l'axe.
event logfile (t += DT_OUTPUT)
{
  static FILE * fp = NULL;
  if (!fp && pid() == 0) {
    fp = fopen ("hauteur_axi_amr.dat", "w");
    if (!fp) { perror ("hauteur_axi_amr.dat"); exit (1); }
    fprintf (fp, "# t_s  h_moyenne_m  h_axe_estime_m\n");
  }

  double volume_reduit = 0., haxe = 0.;
  foreach (reduction(+:volume_reduit) reduction(+:haxe)) {
    volume_reduit += f[]*dv();
    if (y < Delta)
      haxe += f[]*Delta;
  }
  const double hmoy = 2.*volume_reduit/sq(R);
  if (pid() == 0) {
    fprintf (fp, "%.10g %.10g %.10g\n", t,
             hmoy, haxe);
    fflush (fp);
  }
}

// foreach() parcourt les feuilles physiques actives, sans cellules masquees.
event count_cells (t = 0; t += 0.01)
{
  static FILE * fp = NULL;
  if (!fp) {
    fp = fopen ("cellules.dat", "w");
    if (!fp) { perror ("cellules.dat"); exit (1); }
    fprintf (fp, "# t_s cellules_actives\n");
  }
  long ncells = 0;
  foreach (reduction(+:ncells))
    ncells++;
  fprintf (fp, "%.10g %ld\n", t, ncells);
  fflush (fp);
}

// Grille physique fixe aux centres de 1024 x 64 cases, sans points sur les bords.
// interpolate() utilise une interpolation bilineaire et les conditions aux limites.
event export_velocity (t = T_COMPARE)
{
  boundary ((scalar *){u});
  FILE * fp = fopen ("vitesse_t015.dat", "w");
  if (!fp) { perror ("vitesse_t015.dat"); exit (1); }
  fprintf (fp, "# t_s=%.17g\n# x_m y_m ux_m_s\n", t);
  for (int ix = 0; ix < SAMPLE_NX; ix++)
    for (int iy = 0; iy < SAMPLE_NY; iy++) {
      double xp = (ix + 0.5)*HDOM/SAMPLE_NX;
      double yp = (iy + 0.5)*R/SAMPLE_NY;
      double ux = interpolate (u.x, xp, yp);
      if (!isfinite(ux) || ux == nodata) {
        fprintf (stderr, "Interpolation invalide a x=%g, y=%g\n", xp, yp);
        exit (1);
      }
      fprintf (fp, "%.17g %.17g %.17g\n", xp, yp, ux);
    }
  fclose (fp);
}

// Export optionnel active par generer_video_amr.sh.
#if JURIN_EXPORT_VIDEO
#include "jurin_amr_video.h"
#endif

event stop (t = T_END)
{
  return 1;
}

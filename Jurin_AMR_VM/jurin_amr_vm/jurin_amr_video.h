// Export des cellules reelles et des facettes VOF ; calcul en serie.
#ifndef DT_VIDEO
#define DT_VIDEO 0.005
#endif

static double video_last_time = -1.;

void export_video_frame ()
{
  static int frame = 0;
  static FILE * index = NULL;
  char name[80];
  if (!index) {
    FILE * metadata = fopen ("parametres_amr.json", "w");
    if (!metadata) { perror ("parametres_amr.json"); exit (1); }
    fprintf (metadata,
             "{\n  \"radius_m\": %.17g, \"length_m\": %.17g,\n"
             "  \"minlevel\": %d, \"maxlevel\": %d,\n"
             "  \"end_time_s\": %.17g, \"frame_interval_s\": %.17g,\n"
             "  \"height_interval_s\": %.17g, \"f_tolerance\": %.17g,\n"
             "  \"jurin_height_m\": %.17g\n}\n",
             R, HDOM, MINLEVEL, MAXLEVEL, T_END, DT_VIDEO, DT_OUTPUT, F_TOL,
             2.*SIGMA*cos(THETA)/((RHO_L - RHO_G)*G*R));
    fclose (metadata);
    index = fopen ("video_frames.dat", "w");
    if (!index) { perror ("video_frames.dat"); exit (1); }
    fprintf (index, "# frame t_s h_moyenne_m\n");
  }

  sprintf (name, "cells_%04d.dat", frame);
  FILE * fp = fopen (name, "w");
  if (!fp) { perror (name); exit (1); }
  fprintf (fp, "# t=%.17g\n# x_m y_m Delta_m f level\n", t);
  double volume = 0.;
  foreach (serial) {
    fprintf (fp, "%.17g %.17g %.17g %.17g %d\n", x, y, Delta, f[], level);
    volume += f[]*dv();
  }
  fclose (fp);

  sprintf (name, "interface_%04d.dat", frame);
  fp = fopen (name, "w");
  if (!fp) { perror (name); exit (1); }
  output_facets (f, fp);
  fclose (fp);
  fprintf (index, "%d %.17g %.17g\n", frame, t, 2.*volume/sq(R));
  fflush (index);
  if (frame % 20 == 0 || fabs(t - T_END) < 1e-9)
    fprintf (stderr, "Video : t=%.4g / %.4g s (%.0f %%)\n", t, T_END, 100.*t/T_END);
  frame++;
  video_last_time = t;
}

event video_frames (t = 0; t += DT_VIDEO; t <= T_END)
{
  export_video_frame();
}

// L'evenement de fin garantit une image a T_END malgre les arrondis temporels.
event video_last_frame (t = end)
{
  if (fabs(t - video_last_time) > 1e-9)
    export_video_frame();
}

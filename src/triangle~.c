// Copied (by hand) from the cyclone library. There may be typos.

#include "m_pd.h"
#include <string.h>

#define TRIANGLE_DEFPEAK 0.5
#define TRIANGLE_DEFLO -1.0
#define TRIANGLE_DEFHI 1.0

typedef struct _triangle {
  t_object x_obj;
  float x_low;
  float x_range;
  float x_high;
  t_float x_f; // why t_float here and float above?
  t_inlet *x_peaklet;
  t_outlet *x_outlet;
} t_triangle;

static t_class *triangle_class = NULL;

static void triangle_lo(t_triangle *x, t_floatarg f)
{
  x->x_low = f;
  x->x_range = x->x_high - f;
}

static void triangle_hi(t_triangle *x, t_floatarg f)
{
  x->x_high = f;
  x->x_range = f - x->x_low;
}

static t_int *triangle_perform(t_int *w)
{
  t_triangle *x = (t_triangle *)(w[1]);
  int nblock = (int)(w[2]);
  t_float *in1 = (t_float *)(w[3]);
  t_float *in2 = (t_float *)(w[4]);
  t_float *out = (t_float *)(w[5]);
  float low = x->x_low;
  float range = x->x_range;

  while (nblock --) {
    float ph = *in1++;
    float peakph = *in2++;

    if (ph < 0.0) {
      ph -= (int)ph - 1.0; // negative samples get converted to 1
    } else if (ph > 1.0) {
      ph -= (int)ph; // samples > 1 get converted to 0
    }

    if (peakph < 0.0) {
      peakph = 0.0;
    } else if (peakph > 1.0) {
      peakph = 1.0;
    }

    if (ph < peakph) {
      ph /= peakph;
    } else if (peakph < 1.0) {
      ph = (1.0 - ph) / (1.0 - peakph);
    } else {
      ph = 0.0;
    }

    *out++ = low + ph * range;
  }

  return (w + 6);
}

static void triangle_dsp(t_triangle *x, t_signal **sp)
{
  dsp_add(triangle_perform, 5, x, sp[0]->s_length,
          sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec);
}

static void *triangle_new(t_symbol *s, int argc, t_atom *argv)
{
  t_triangle *x = (t_triangle *)pd_new(triangle_class);

  t_float tripeak = TRIANGLE_DEFPEAK;
  t_float trilo = x->x_low = TRIANGLE_DEFLO;
  t_float trihi = x->x_high = TRIANGLE_DEFHI;

  int argnum = 0;
  while(argc > 0) {
    if (argv->a_type == A_FLOAT) { // if current arg is a number
      t_float argval = atom_getfloatarg(0, argc, argv);
      switch(argnum) {
        case 0:
          tripeak = argval;
          break;
        default:
          break;
      }
      argnum++;
      argc--;
      argv++;
    } else if (argv->a_type == A_SYMBOL) {
      t_symbol *curarg = atom_getsymbolarg(0, argc, argv);
      if (strcmp(curarg->s_name, "@lo") == 0) {
        if (argc >= 2) {
          trilo = atom_getfloatarg(1, argc, argv);
          argc -= 2;
          argv += 2;
        } else {
          goto errstate;
        }
      } else if (strcmp(curarg->s_name, "@hi") == 0) {
        if (argc >= 2) {
          trihi = atom_getfloatarg(1, argc, argv);
          argc -= 2;
          argv += 2;
        } else {
          goto errstate;
        }
      } else {
        goto errstate;
      }
    } else {
      goto errstate;
    }
  }

  triangle_lo(x, trilo);
  triangle_hi(x, trihi);

  x->x_peaklet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_peaklet, tripeak);
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  return (void *)x;
errstate:
  pd_error(x, "triangle~: improper args");
  return NULL;
}

void *triangle_free(t_triangle *x)
{
  inlet_free(x->x_peaklet);
  outlet_free(x->x_outlet);

  return (void *)x;
}

void triangle_tilde_setup(void)
{
  triangle_class = class_new(gensym("triangle~"),
                             (t_newmethod)triangle_new,
                             (t_method)triangle_free,
                             sizeof(t_triangle),
                             CLASS_DEFAULT,
                             A_GIMME, 0);

  class_addmethod(triangle_class, (t_method)triangle_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(triangle_class, t_triangle, x_f);
  class_addmethod(triangle_class, (t_method)triangle_lo,
                  gensym("lo"), A_DEFFLOAT, 0);
  class_addmethod(triangle_class, (t_method)triangle_hi,
                  gensym("hi"), A_DEFFLOAT, 0);
}


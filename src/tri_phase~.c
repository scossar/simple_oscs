#include "m_pd.h"

static t_class *tri_phase_class = NULL;

/*
 * See `simple_phasor.c` for detailed notes
 */
#define UNITBIT32 1572864.  /* 3*2^19; bit 32 has place value 1 */

#if BYTE_ORDER == LITTLE_ENDIAN
# define HIOFFSET 1
# define LOWOFFSET 0
#else
# define HIOFFSET 0    /* word offset to find MSB */
# define LOWOFFSET 1    /* word offset to find LSB */
#endif

union tabfudge
{
    double tf_d;
    int32_t tf_i[2];
};

typedef struct _tri_phase
{
  t_object x_obj;
  double x_phase;
  t_float x_conv;
  t_float x_f;

  t_float x_low;
  t_float x_hi;
  t_float x_peak;
  t_float x_range;

  t_inlet *in_1; // frequency
  t_inlet *in_2; // peak
  t_inlet *in_3; // phase
} t_tri_phase;


static t_int *tri_phase_perform(t_int *w)
{
  t_tri_phase *x = (t_tri_phase *)(w[1]);
  t_sample *in1 = (t_float *)(w[2]);
  t_sample *in2 = (t_float *)(w[3]);
  t_sample *out = (t_float *)(w[4]);
  int n = (int)(w[5]);
  double dphase = x->x_phase + (double)UNITBIT32;
  union tabfudge tf;
  int normhipart;
  t_float conv = x->x_conv;

  // hardcoded for now
  float low = x->x_low;
  float range = x->x_range;

  tf.tf_d = UNITBIT32;
  normhipart = tf.tf_i[HIOFFSET];
  tf.tf_d = dphase;

  while (n--)
  {
    float peak = *in2++;

    if (peak < 0.0) {
      peak = 0.0;
    } else if (peak > 1.0) {
      peak = 1.0;
    }

    tf.tf_i[HIOFFSET] = normhipart;
    dphase += *in1++ * conv;

    float ph = tf.tf_d - UNITBIT32;

    if (ph < peak) {
      ph /= peak;
    } else if (peak < 1.0) {
      ph = (1.0 - ph) / (1.0 - peak);
    } else {
      ph = 0.0;
    }

    *out++ = low + ph * range;
    tf.tf_d = dphase;
  }

  tf.tf_i[HIOFFSET] = normhipart;
  x->x_phase = tf.tf_d - UNITBIT32;
  return (w+6);
}

static void tri_phase_dsp(t_tri_phase *x, t_signal **sp)
{
  x->x_conv = 1.0 / sp[0]->s_sr;
  dsp_add(tri_phase_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, (t_int)sp[0]->s_length);
}

static void tri_phase_ft1(t_tri_phase *x, t_float f)
{
  x->x_phase = (double)f;
}

static void tri_phase_low(t_tri_phase *x, t_floatarg f)
{
  x->x_low = f;
  x->x_range = x->x_hi - f;
}

static void tri_phase_hi(t_tri_phase *x, t_floatarg f)
{
  x->x_hi = f;
  x->x_range = f - x->x_low;
}

static void *tri_phase_new(t_floatarg f)
{
  t_tri_phase *x = (t_tri_phase *)pd_new(tri_phase_class);
  x->x_f = f;
  x->x_phase = 0;
  x->x_conv = 0;
  x->x_low = -1.0;
  x->x_hi = 1.0;
  x->x_peak = 0.5;

  tri_phase_low(x, x->x_low);
  tri_phase_hi(x, x->x_hi);

  x->in_1 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->in_1, x->x_f);

  x->in_2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->in_2, x->x_peak);

  x->in_3 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("ft1"));

  outlet_new(&x->x_obj, gensym("signal"));

  return (void *)x;
}

void tri_phase_tilde_setup(void)
{
  tri_phase_class = class_new(gensym("tri_phase~"), (t_newmethod)tri_phase_new, 0,
                                  sizeof(t_tri_phase), 0, A_DEFFLOAT, 0);
  CLASS_MAINSIGNALIN(tri_phase_class, t_tri_phase, x_f);
  class_addmethod(tri_phase_class, (t_method)tri_phase_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(tri_phase_class, (t_method)tri_phase_ft1,
                  gensym("ft1"), A_FLOAT, 0);
}

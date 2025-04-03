#include "m_pd.h"
#include <math.h>

static t_class *tri_phase_class = NULL;

/*
 * See `simple_phasor.c` for detailed notes
 */
#define UNITBIT32 1572864.  /* 3*2^19; bit 32 has place value 1 */

// I suspect this is already set via the Pure Data environment.
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
  t_float x_softness;

  t_float fold_threshold;

  t_inlet *in_2; // peak
  t_inlet *in_3; // fold_threshold
  t_inlet *in_4; // fold softness
  t_inlet *in_5; // phase
} t_tri_phase;

static float tri_phase_fold(float sample, float threshold, float softness)
{
  float transition_region = softness * threshold;
  if (transition_region > 0.0f) {
    while (1) {
      if (sample > threshold) {
        float overshoot = sample - threshold;
        if (overshoot <= transition_region) {
          // apply cubic curve in transition region
          float t = overshoot / transition_region;
          float fold_amount = t*t * (3.0f - 2.0f*t); // cubic hermite
          sample = threshold - fold_amount * overshoot;
          break;
        } else {
          // regular fold, but account for transition region
          sample = 2.0f * threshold - sample;
        }
      } else if (sample < -threshold)
      {
        float overshoot = -threshold - sample;
        if (overshoot <= transition_region) {
          float t = overshoot / transition_region;
          float fold_amount = t*t * (3.0f - 2.0f*t); // cubic hermite
          sample = -threshold + fold_amount * overshoot;
          break;
        } else {
          // regular fold, but account for transition region
          sample = -2.0f * threshold - sample;
        }
      } else {
        break; // no folding needed
      }

    }
  } else {
    // regular hard folding
    while (sample > threshold) {
      sample = 2.0f * threshold - sample;
    }
    while (sample < -threshold) {
      sample = -2.0f * threshold - sample;
    }
  }

  return sample;
}


static t_int *tri_phase_perform(t_int *w)
{
  t_tri_phase *x = (t_tri_phase *)(w[1]);
  t_sample *in1 = (t_float *)(w[2]); // frequency input
  t_sample *in2 = (t_float *)(w[3]); // peak input
  t_sample *in3 = (t_float *)(w[4]); // fold threshold input
  t_sample *out = (t_float *)(w[5]); // output
  int n = (int)(w[6]);

  double dphase = x->x_phase + (double)UNITBIT32;
  union tabfudge tf;
  int normhipart;
  t_float conv = x->x_conv;
  // hardcoded for now
  float low = x->x_low;
  float range = x->x_range;

  float softness = x->x_softness;

  tf.tf_d = UNITBIT32;
  normhipart = tf.tf_i[HIOFFSET];
  tf.tf_d = dphase;

  while (n--)
  {
    float peak = *in2++;
    peak = (peak < 0.0f) ? 0.0f : (peak > 1.0f) ? 1.0f : peak;

    float threshold = *in3++;
    threshold = (threshold < 0.0f) ? 0.0f : threshold;

    // update phase
    tf.tf_i[HIOFFSET] = normhipart;
    dphase += *in1++ * conv;
    float ph = tf.tf_d - UNITBIT32;

    // generate triangle wave with variable peak
    float tri_value;
    if (ph < peak) {
      tri_value = (peak > 0.0f) ? ph / peak : 0.0f;
    } else if (peak < 1.0f) {
      tri_value = (1.0f - peak > 0.0f) ? (1.0f - ph) / (1.0f - peak) : 0.0f;
    } else {
      tri_value = 0.0f;
    }

    // scale to output range
    float s = low + tri_value * range;

    // apply wave folding
    s = tri_phase_fold(s, threshold, softness);

    *out++ = s;
    tf.tf_d = dphase;
  }

  tf.tf_i[HIOFFSET] = normhipart;
  x->x_phase = tf.tf_d - UNITBIT32;
  return (w+7);
}

static void tri_phase_dsp(t_tri_phase *x, t_signal **sp)
{
  x->x_conv = 1.0 / sp[0]->s_sr;
  dsp_add(tri_phase_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, (t_int)sp[0]->s_length);
}

static void tri_phase_ft1(t_tri_phase *x, t_float f)
{
  x->x_phase = (double)f;
}

static void tri_phase_softness(t_tri_phase *x, t_float f)
{
  x->x_softness = f;
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

static void tri_phase_free(t_tri_phase *x)
{
  inlet_free(x->in_2);
  inlet_free(x->in_3);
  inlet_free(x->in_4);
  inlet_free(x->in_5);
}

static void *tri_phase_new(t_floatarg f)
{
  t_tri_phase *x = (t_tri_phase *)pd_new(tri_phase_class);
  x->x_f = f;
  x->x_phase = 0;
  x->x_conv = 0;

  // TODO: make configurable
  x->x_low = -1.0;
  x->x_hi = 1.0;

  x->fold_threshold = 0.5; // default
  x->x_peak = 0.5; // default
  x->x_softness = 0.5;

  tri_phase_low(x, x->x_low);
  tri_phase_hi(x, x->x_hi);

  // where the peak is in the phase (0, 1)
  x->in_2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->in_2, x->x_peak);

  // threshold for folding (0, 1)
  x->in_3 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->in_3, x->fold_threshold);

  // folding softness
  x->in_4 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("softness"));

  x->in_5 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("ft1"));

  outlet_new(&x->x_obj, gensym("signal"));

  return (void *)x;
}

void tri_phase_tilde_setup(void)
{
  tri_phase_class = class_new(gensym("tri_phase~"),
                              (t_newmethod)tri_phase_new,
                              (t_method)tri_phase_free,
                              sizeof(t_tri_phase),
                              CLASS_DEFAULT,
                              A_DEFFLOAT, 0);

  CLASS_MAINSIGNALIN(tri_phase_class, t_tri_phase, x_f);
  class_addmethod(tri_phase_class, (t_method)tri_phase_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(tri_phase_class, (t_method)tri_phase_ft1,
                  gensym("ft1"), A_FLOAT, 0);
  class_addmethod(tri_phase_class, (t_method)tri_phase_softness,
                  gensym("softness"), A_FLOAT, 0);
}

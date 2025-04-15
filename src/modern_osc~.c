// updated version of simple_osc~.c

#include "m_pd.h"
#include <math.h>

// #define WAVETABLE_SIZE 16384 // 2^14
#define WAVETABLE_SIZE 4096 // 2^12 might be good enough

static t_class *modern_osc_class = NULL;
static float *cos_table = NULL; // shared wavetable
static int table_reference_count = 0; // track how many instances exist

typedef struct _modern_osc {
  t_object x_obj;
  double x_phase;
  t_float x_conv;
  t_inlet *x_freq_inlet;
  t_outlet *x_outlet;
  t_float x_f;
} t_modern_osc;

static void wavetable_init(void)
{
  if (cos_table == NULL) {
    cos_table = (float *)getbytes(sizeof(float) * (WAVETABLE_SIZE ));
    if (cos_table) {
      for (int i = 0; i < WAVETABLE_SIZE; i++) {
        cos_table[i] = cosf((i * 2.0f * (float)M_PI) / (float)WAVETABLE_SIZE);
      }
      post("modern_osc~: initialized cosine table of size %d", WAVETABLE_SIZE);
    } else {
      post("modern_osc~ error: failed to allocate memory for cosine table");
    }
  }
  table_reference_count++;
}

static void wavetable_free(void)
{
  table_reference_count--;
  if (table_reference_count <= 0 && cos_table != NULL) {
    freebytes(cos_table, sizeof(float) * (WAVETABLE_SIZE));
    cos_table = NULL;
    post("modern_osc~: freed cosine table");
    table_reference_count = 0; // just to be safe
  }
}

static t_int *modern_osc_perform(t_int *w)
{
  t_modern_osc *x = (t_modern_osc *)(w[1]);
  t_sample *in = (t_sample *)(w[2]); // fix the type
  t_sample *out = (t_sample *)(w[3]); // fix the type
  int n = (int)(w[4]);

  float *tab = cos_table;
  t_float conv = x->x_conv;
  double phase = x->x_phase;

  if (!tab) return (w+5);

  while (n--) {
    double curphase = phase;
    phase += *in++ * conv;
    unsigned int idx = (unsigned int)curphase;
    t_sample frac = (t_sample)(curphase - idx);

    idx &= (WAVETABLE_SIZE - 1);

    t_sample f1 = tab[idx];
    t_sample f2 = tab[(idx + 1) & (WAVETABLE_SIZE - 1)];
    *out++ = f1 + frac * (f2 - f1);
  }

  while (phase >= WAVETABLE_SIZE) phase -= WAVETABLE_SIZE;
  while (phase < 0) phase += WAVETABLE_SIZE;
  x->x_phase = phase;

  return (w + 5);
}

static void modern_osc_dsp(t_modern_osc *x, t_signal **sp)
{
  // calculate the conversion factor for this sample rate
  x->x_conv = (float)WAVETABLE_SIZE / sp[0]->s_sr;

  dsp_add(modern_osc_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
}

static void *modern_osc_new(t_floatarg f)
{
  t_modern_osc *x = (t_modern_osc *)pd_new(modern_osc_class);

  x->x_phase = (double)0.0;
  x->x_f = f > 0 ? (t_float)f : (t_float)220.0;

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f);
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  wavetable_init();

  return (void *)x;
}

static void modern_osc_free(t_modern_osc *x)
{
  if (x->x_freq_inlet) {
    inlet_free(x->x_freq_inlet);
  }

  if (x->x_outlet) {
    outlet_free(x->x_outlet);
  }

  // decrease reference count and possibly free wavetable
  wavetable_free();
}

void modern_osc_tilde_setup(void)
{
  modern_osc_class = class_new(gensym("modern_osc~"),
                               (t_newmethod)modern_osc_new,
                               (t_method)modern_osc_free,
                               sizeof(t_modern_osc),
                               CLASS_DEFAULT,
                               A_DEFFLOAT, 0);

  class_addmethod(modern_osc_class, (t_method)modern_osc_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(modern_osc_class, t_modern_osc, x_f);
}



#include "m_pd.h"
#include <math.h>

// #define WAVETABLE_SIZE 16384
#define WAVETABLE_SIZE 65536

static t_class *fold_osc_class = NULL;
static t_float *cos_table = NULL; // shared wavetable
static int table_reference_count = 0; // track how many instances exist

typedef struct _fold_osc {
  t_object x_obj;
  double x_phase;
  t_float x_conv;
  t_inlet *x_freq_inlet;
  t_inlet *x_fold_inlet;
  t_outlet *x_outlet;
  t_float x_f;
  t_float x_threshold; // fold threshold value
} t_fold_osc;

static void wavetable_init(void)
{
  if (cos_table == NULL) {
    cos_table = (t_float *)getbytes(sizeof(t_float) * (WAVETABLE_SIZE + 1));
    if (cos_table) {
      for (int i = 0; i <= WAVETABLE_SIZE; i++) {
        cos_table[i] = cosf((i * 2.0f * M_PI) / WAVETABLE_SIZE);
      }
      post("fold_osc~: initialized cosine table of size %d", WAVETABLE_SIZE);
    } else {
      post("fold_osc~ error: failed to allocate memory for cosine table");
    }
  }
  table_reference_count++;
}

static void wavetable_free(void)
{
  table_reference_count--;
  if (table_reference_count <= 0 && cos_table != NULL) {
    freebytes(cos_table, sizeof(t_float) * (WAVETABLE_SIZE + 1));
    cos_table = NULL;
    post("simple_osc~: freed cosine table");
    table_reference_count = 0; // just to be safe
  }
}

static float cubicInterpolate(float y0, float y1, float y2, float y3, float mu) {
    float a0, a1, a2, a3, mu2;
    
    mu2 = mu * mu;
    a0 = y3 - y2 - y0 + y1;
    a1 = y0 - y1 - a0;
    a2 = y2 - y0;
    a3 = y1;
    
    return a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3;
}

static t_int *fold_osc_perform(t_int *w)
{
  t_fold_osc *x = (t_fold_osc *)(w[1]);
  t_float *in1 = (t_float *)(w[2]);
  t_float *in2 = (t_float *)(w[3]);
  t_float *out = (t_float *)(w[4]);
  int n = (int)(w[5]);

  double dphase = x->x_phase;
  double conv = x->x_conv;
  // t_float threshold = x->x_threshold;

  if (!cos_table) return (w+5);

  while (n--) {
    t_float freq = *in1++;
    t_float current_threshold = *in2++;
    // generate two samples per output sample (2x oversampling)

    t_float sample = 0.0f;
    for (int i = 0; i < 2; i++) {
      int index = ((int)dphase) & (WAVETABLE_SIZE - 1);
      t_float frac = dphase - index; // the fractional part after getting the
      // index int
      t_float y0 = cos_table[(index - 1) & (WAVETABLE_SIZE - 1)]; // in case
      // index - 1 is out of range
      t_float y1 = cos_table[index];
      t_float y2 = cos_table[(index + 1) & (WAVETABLE_SIZE - 1)];
      t_float y3 = cos_table[(index + 2) & (WAVETABLE_SIZE - 1)];

      t_float oscillator_out = cubicInterpolate(y0, y1, y2, y3, frac);
      if (oscillator_out > current_threshold) {
        oscillator_out = 2.0f * current_threshold - oscillator_out;
      } else if (oscillator_out < -current_threshold) {
        oscillator_out = -2.0f * current_threshold - oscillator_out;
      }

      // accumulate for averaging
      sample += oscillator_out * 0.5f;

      // advance phase at half the increment for oversampling
      dphase += (freq * conv) * 0.5f;
      while (dphase >= WAVETABLE_SIZE) dphase -= WAVETABLE_SIZE;
      while (dphase < 0) dphase += WAVETABLE_SIZE;
    }

    *out++ = sample;
  }

  x->x_phase = dphase;
  return (w + 6);
}

static void fold_osc_dsp(t_fold_osc *x, t_signal **sp)
{
  // calculate the conversion factor for this sample rate
  x->x_conv = WAVETABLE_SIZE / sp[0]->s_sr;

  dsp_add(fold_osc_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_length);
}

static void *fold_osc_new(t_floatarg f)
{
  t_fold_osc *x = (t_fold_osc *)pd_new(fold_osc_class);

  // initialize phase and frequency
  x->x_phase = 0;
  x->x_f = f > 0 ? f : 440;
  x->x_threshold = 0.5f;

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f); // sets inlet initial value from
  // arg

  x->x_fold_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_fold_inlet, x->x_threshold);

  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  wavetable_init();

  return (void *)x;
}

static void fold_osc_free(t_fold_osc *x)
{
  inlet_free(x->x_freq_inlet);
  inlet_free(x->x_fold_inlet);
  outlet_free(x->x_outlet);

  // decrease reference count and possibly free wavetable
  wavetable_free();
}

void fold_osc_tilde_setup(void)
{
  fold_osc_class = class_new(gensym("fold_osc~"),
                               (t_newmethod)fold_osc_new,
                               (t_method)fold_osc_free,
                               sizeof(t_fold_osc),
                               CLASS_DEFAULT,
                               A_DEFFLOAT, 0);

  class_addmethod(fold_osc_class, (t_method)fold_osc_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(fold_osc_class, t_fold_osc, x_f);
}

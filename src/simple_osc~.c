// reference: pd_cyclone_triangle_external.md (line: ~1020)
// for ideas about converting this to a wave folding oscillator, see
// pure_data_wave_folding_oscillator.md
//
// for a general overview of the code, (mostly mine with some confirmation from
// Claude), see `understanding_pure_data_wave_table_oscillator.md`

#include "m_pd.h"
#include <math.h>

#define WAVETABLE_SIZE 16384

static t_class *simple_osc_class = NULL;
static t_float *cos_table = NULL; // shared wavetable
static int table_reference_count = 0; // track how many instances exist

typedef struct _simple_osc {
  t_object x_obj;
  double x_phase;
  t_float x_conv;
  t_inlet *x_freq_inlet;
  t_outlet *x_outlet;
  t_float x_f;
} t_simple_osc;

static void wavetable_init(void)
{
  if (cos_table == NULL) {
    cos_table = (t_float *)getbytes(sizeof(t_float) * (WAVETABLE_SIZE + 1));
    if (cos_table) {
      for (int i = 0; i <= WAVETABLE_SIZE; i++) {
        cos_table[i] = cosf((i * 2.0f * M_PI) / WAVETABLE_SIZE);
      }
      post("simple_osc~: initialized cosine table of size %d", WAVETABLE_SIZE);
    } else {
      post("simple_osc~ error: failed to allocate memory for cosine table");
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

static t_int *simple_osc_perform(t_int *w)
{
  t_simple_osc *x = (t_simple_osc *)(w[1]);
  t_float *in = (t_float *)(w[2]);
  t_float *out = (t_float *)(w[3]);
  int n = (int)(w[4]);

  double dphase = x->x_phase;
  double conv = x->x_conv;

  if (!cos_table) return (w+5);

  while (n--) {
    t_float freq = *in++;
    int index = ((int)dphase) & (WAVETABLE_SIZE-1);
    t_float frac = dphase - index;

    // linear interpolation between table points
    *out++ = cos_table[index] + frac * (cos_table[index + 1] - cos_table[index]);

    // advance phase based on frequency
    dphase += freq * conv;
    while (dphase >= WAVETABLE_SIZE) dphase -= WAVETABLE_SIZE;
    while (dphase <0) dphase += WAVETABLE_SIZE;
  }

  x->x_phase = dphase;
  return (w + 5);
}

static void simple_osc_dsp(t_simple_osc *x, t_signal **sp)
{
  // calculate the conversion factor for this sample rate
  x->x_conv = WAVETABLE_SIZE / sp[0]->s_sr;

  dsp_add(simple_osc_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
}

static void *simple_osc_new(t_floatarg f)
{
  t_simple_osc *x = (t_simple_osc *)pd_new(simple_osc_class);

  // initialize phase and frequency
  x->x_phase = 0;
  x->x_f = f > 0 ? f : 440;

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f); // sets inlet initial value from
  // arg
  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  wavetable_init();

  return (void *)x;
}

static void simple_osc_free(t_simple_osc *x)
{
  inlet_free(x->x_freq_inlet);
  outlet_free(x->x_outlet);

  // decrease reference count and possibly free wavetable
  wavetable_free();
}

void simple_osc_tilde_setup(void)
{
  simple_osc_class = class_new(gensym("simple_osc~"),
                               (t_newmethod)simple_osc_new,
                               (t_method)simple_osc_free,
                               sizeof(t_simple_osc),
                               CLASS_DEFAULT,
                               A_DEFFLOAT, 0);

  class_addmethod(simple_osc_class, (t_method)simple_osc_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(simple_osc_class, t_simple_osc, x_f);
}



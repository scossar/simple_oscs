// reference pure_data_osc_perform.md

#include "m_pd.h"
#include <math.h>

// I'm not sure the table needs to be so big. It does need to be a power of 2
// though
#define WAVETABLE_SIZE 16384 // 2^14

// 3*2^19
// When added to a fractional value between 0 and 1, and assigned as a double, the integer part of the
// result is stored in the high word of the double. This will be used to get the
// index of the cosine table.
#define UNITBIT32 1572864.0

#if BYTE_ORDER == LITTLE_ENDIAN
# define HIOFFSET 1
# define LOWOFFSET 0
#else
# define HIOFFSET 0    /* word offset to find MSB */
# define LOWOFFSET 1    /* word offset to find LSB */
#endif

static t_class *tabfudge_osc_class = NULL;
// using `float` intentionally here, see pd_floattype.md
static float *cos_table = NULL; 
static int table_reference_count = 0; // tracks shared instances of cos_table

union tabfudge {
  double tf_d;
  int32_t tf_i[2];
};

typedef struct _tabfudge_osc {
  t_object x_obj;
  double x_phase;
  t_float x_conv;
  t_inlet *x_freq_inlet;
  t_outlet *x_outlet;
  t_float x_f;
} t_tabfudge_osc;


static void wavetable_init(void)
{
  if (cos_table == NULL) {
    // NOTE: I'm unsure/undecided about WAVETABLE_SIZE + 1; it makes the perform
    // function a little simpler
    cos_table = (float *)getbytes(sizeof(float) * (WAVETABLE_SIZE + 1));
    if (cos_table) {
      for (int i = 0; i <= WAVETABLE_SIZE; i++) { // <= account for the +1
        cos_table[i] = cosf((i * 2.0f * (float)M_PI) / (float)WAVETABLE_SIZE);
      }
      post("tabfudge_osc~: initialized cosine table of size %d", WAVETABLE_SIZE);
    } else {
      post("tabfudge_osc~ error: failed to allocate memory for cosine table");
    }
  }
  table_reference_count++;
}

static void wavetable_free(void)
{
  table_reference_count--;
  if (table_reference_count <= 0 && cos_table != NULL) {
    // NOTE: if + 1 is removed from the init function, remove it here as well
    freebytes(cos_table, sizeof(float) * (WAVETABLE_SIZE + 1));
    cos_table = NULL;
    post("tabfudge_osc~: freed cosine table");
    table_reference_count = 0; // just to be safe
  }
}

static t_int *tabfudge_osc_perform(t_int *w)
{
  t_tabfudge_osc *x = (t_tabfudge_osc *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *out1 = (t_sample *)(w[3]);
  int n = (int)(w[4]);

  float *tab = cos_table;
  float *addr;
  t_sample f1, f2, frac;
  // when assigned to tf.tf_d:
  // - the high word of this double will be the wavetable index
  // - subtracting UNITBIT32 will give the fractional part of the index
  double dphase = x->x_phase + UNITBIT32;
  int normhipart;
  union tabfudge tf;
  float conv = x->x_conv;

  if (!tab) return (w+5);

  tf.tf_d = UNITBIT32; // initialize union 
  normhipart = tf.tf_i[HIOFFSET]; // save the high word of (double)UNITBIT32

  while (n--) {
    tf.tf_d = dphase; // see dphase comment
    // update dphase with freq_input * sample rate conversion
    dphase += *in1 * conv;
    // pointer to wavetable + index % wavetable (using bit mask)
    // performs index extraction and modulo operation in 1 step
    addr = tab + (tf.tf_i[HIOFFSET] & (WAVETABLE_SIZE - 1));
    // reset hi word to isolate fractional part
    tf.tf_i[HIOFFSET] = normhipart;
    // extract just the fractional part by subtracting UNITBIT32
    frac = tf.tf_d - UNITBIT32;
    f1 = addr[0];
    // the +1 in wavetable_init allows for this
    // the last wavetable memory slot should hold an equal value to the first
    // slot
    f2 = addr[1];
    // interpolation; unsure why amplitudes don't need to be scaled
    *out1++ = f1 + frac * (f2 - f1);
  }

  // oh no... 
  // the code below is a more effecient version of
  // x->x_phase = fmod(dphase - UNITBIT32, WAVETABLE_SIZE);
  //
  // establish a reference value that represents one cycle in the bit-encoded
  // number space
  // this is a more efficient version of
  tf.tf_d = UNITBIT32 * WAVETABLE_SIZE;
  // the high word of the value is stored in normhipart
  normhipart = tf.tf_i[HIOFFSET];
  // set tf.tf_d again, this time as current dphase plus one complete cycle
  // minus UNITBIT32
  tf.tf_d = dphase + (UNITBIT32 * WAVETABLE_SIZE - UNITBIT32);
  // reset HIOFFSET to the normhipart of the reference value
  // effectively performing a modulo operation on the phase
  tf.tf_i[HIOFFSET] = normhipart;
  // subtract a full table size from the value
  // results in a value in the range [0, WAVETABLE_SIZE]
  x->x_phase = tf.tf_d - UNITBIT32 * WAVETABLE_SIZE;

#if 0 // human readable equivalent
while (dphase >= WAVETABLE_SIZE) dphase -= WAVETABLE_SIZE;
while (dphase < 0) dphase += WAVETABLE_SIZE;
#endif

  return (w+5);
}

static void tabfudge_osc_dsp(t_tabfudge_osc *x, t_signal **sp)
{
  x->x_conv = (float)WAVETABLE_SIZE / sp[0]->s_sr;
  dsp_add(tabfudge_osc_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
}

static void tabfudge_osc_free(t_tabfudge_osc *x)
{
  if (x->x_freq_inlet) {
    inlet_free(x->x_freq_inlet);
  }

  if (x->x_outlet) {
    outlet_free(x->x_outlet);
  }

  wavetable_free();
}

static void *tabfudge_osc_new(t_floatarg f)
{
  t_tabfudge_osc *x = (t_tabfudge_osc *)pd_new(tabfudge_osc_class);

  x->x_f = f > 0 ? (t_float)f : (t_float)220.0;
  x->x_phase = (double)0.0;

  x->x_freq_inlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
  pd_float((t_pd *)x->x_freq_inlet, x->x_f);

  x->x_outlet = outlet_new(&x->x_obj, &s_signal);

  wavetable_init();

  return (void *)x;
}

void tabfudge_osc_tilde_setup(void)
{
  tabfudge_osc_class = class_new(gensym("tabfudge_osc~"),
                                 (t_newmethod)tabfudge_osc_new,
                                 (t_method)tabfudge_osc_free,
                                 sizeof(t_tabfudge_osc),
                                 CLASS_DEFAULT,
                                 A_DEFFLOAT, 0);

  class_addmethod(tabfudge_osc_class, (t_method)tabfudge_osc_dsp, gensym("dsp"), A_CANT, 0);
  CLASS_MAINSIGNALIN(tabfudge_osc_class, t_tabfudge_osc, x_f);
}




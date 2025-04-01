#include "m_pd.h"

static t_class *simple_phasor_class = NULL;

/*
 * See `pure_data_phasor_external.md` for notes about what's going on with
 * UNITBIT32 and tabfudge. It's a bit of magic that allows the fractional
 * (in the range (0, 1)) phase values to be wrapped without using a modulo
 * operator.
 */
#define UNITBIT32 1572864.  /* 3*2^19; bit 32 has place value 1 */

// I think this is handled elsewhere, but just in case:
// BYTE_ORDER should be defined by the operation system. Most modern systems are
// little endian (least significant bytes stored first (remember bytes not bits!))
#if BYTE_ORDER == LITTLE_ENDIAN
# define HIOFFSET 1
# define LOWOFFSET 0
#else
# define HIOFFSET 0    /* word offset to find MSB */
# define LOWOFFSET 1    /* word offset to find LSB */
#endif

// Overlayes a double (64 bits) with an array of two tf_i values (each has 32
// bits). 
// tf_i[HIOFFSET] contains the sign bit, the exponent bits, and the first 20
// bits of the mantissa (fractional part).
// tf_i[LOWOFFSET] stores the rest of the mantissa.
// Because it's a union, both members share the same memory location - updating
// tf_i[HIOFFSET] in the while loop modifies the corresponding 32 bits in the 64
// bit tf_d double.
union tabfudge
{
    double tf_d;
    int32_t tf_i[2];
};

typedef struct _simple_phasor
{
  t_object x_obj;
  double x_phase;
  t_float x_conv;
  t_float x_f; // scalar frequency
} t_simple_phasor;

static void *simple_phasor_new(t_floatarg f)
{
  t_simple_phasor *x = (t_simple_phasor *)pd_new(simple_phasor_class);
  x->x_f = f;
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("ft1"));
  x->x_phase = 0;
  x->x_conv = 0;
outlet_new(&x->x_obj, gensym("signal"));

  return (void *)x;
}

static t_int *simple_phasor_perform(t_int *w)
{
  t_simple_phasor *x = (t_simple_phasor *)(w[1]);
  t_sample *in = (t_float *)(w[2]);
  t_sample *out = (t_float *)(w[3]);
  int n = (int)(w[4]);
  double dphase = x->x_phase + (double)UNITBIT32;
  union tabfudge tf;
  int normhipart;
  t_float conv = x->x_conv;

  tf.tf_d = UNITBIT32;
  normhipart = tf.tf_i[HIOFFSET]; // get the HIOFFSET of UNITBIT32. The binary
  // representation of 3*2^19.
  tf.tf_d = dphase; // overwrite tf.tf_d to dphase + UNITBIT32. The high part
  // now matches UNITBIT32 and the low part carries the fractional value from
  // dphase

  while (n--)
  {
    tf.tf_i[HIOFFSET] = normhipart; // keep setting the hioffset to normhipart
    dphase += *in++ * conv; // frequency * conv factor
    *out++ = tf.tf_d - UNITBIT32; // just output the fractional part
    tf.tf_d = dphase; // reset for next iteration
  }

  tf.tf_i[HIOFFSET] = normhipart;
  x->x_phase = tf.tf_d - UNITBIT32;
  return (w+5);
}

static void simple_phasor_dsp(t_simple_phasor *x, t_signal **sp)
{
  x->x_conv = 1./sp[0]->s_sr;
  dsp_add(simple_phasor_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, (t_int)sp[0]->s_length);
}

static void simple_phasor_ft1(t_simple_phasor *x, t_float f)
{
  x->x_phase = (double)f;
}

void simple_phasor_tilde_setup(void)
{
  simple_phasor_class = class_new(gensym("simple_phasor~"), (t_newmethod)simple_phasor_new, 0,
                                  sizeof(t_simple_phasor), 0, A_DEFFLOAT, 0);
  CLASS_MAINSIGNALIN(simple_phasor_class, t_simple_phasor, x_f);
  class_addmethod(simple_phasor_class, (t_method)simple_phasor_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(simple_phasor_class, (t_method)simple_phasor_ft1,
                  gensym("ft1"), A_FLOAT, 0);
}

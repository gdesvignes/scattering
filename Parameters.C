#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "Parameters.h"

int readParameters(param *p, char *paramFile){
    GKeyFile* gkf; /* Notice we declared a pointer */

    char string[32];

    gkf = g_key_file_new();
    gsize length;

    if (!g_key_file_load_from_file(gkf, paramFile, G_KEY_FILE_NONE, NULL)){
        fprintf (stderr, "Could not read config file %s\n", paramFile);
        return EXIT_FAILURE;
    }

    p->IS = g_key_file_get_integer(gkf,"multinest","IS",NULL);
    p->nlive = g_key_file_get_integer(gkf,"multinest","nlive",NULL);
    p->ceff = g_key_file_get_integer(gkf,"multinest","ceff",NULL);
    p->efr = g_key_file_get_double(gkf,"multinest","efr",NULL);
    if (g_key_file_has_key(gkf,"multinest","basename",NULL))
        p->basename = g_key_file_get_string(gkf,"multinest","basename",NULL);

    //p->threshold = g_key_file_get_double(gkf,"config","threshold",NULL);
    //p->have_efac = g_key_file_get_integer(gkf,"config","have_efac",NULL);
    //p->margin_phi0 = g_key_file_get_integer(gkf,"config","margin_phi0",NULL);

    p->sigma = g_key_file_get_double_list(gkf,"params","sigma", &length, NULL);
    p->tau = g_key_file_get_double_list(gkf,"params","tau", &length, NULL);
    p->t0_inf = g_key_file_get_double_list(gkf,"params","t0inf", &length, NULL);
    p->amp = g_key_file_get_double_list(gkf,"params","amp", &length, NULL);
    p->b = g_key_file_get_double_list(gkf,"params","baseline", &length, NULL);
    p->DM = g_key_file_get_double_list(gkf,"params","DM", &length, NULL);
    p->nchan = g_key_file_get_integer(gkf,"params","nchan", NULL);
    p->bscrunch = g_key_file_get_integer(gkf,"params","bscrunch", NULL);

    fprintf (stderr, "Finished reading config file %s\n", paramFile);

    g_key_file_free (gkf);

    return EXIT_SUCCESS;
}

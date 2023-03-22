#include <iostream>
#include "lpmd_calc.cpp"
#include "array_dist_calc.cpp"

static void usage(){
    fprintf(fp, "\n");
    fprintf(fp, "Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)\n");
#if USE_GPL
    fprintf(fp, "License: GNU GPLv3+, due to use of the GNU Scientific Library\n");
#endif
    fprintf(fp, "Version: %s (using htslib %s)\n", bcftools_version(), hts_version());
    fprintf(fp, "\n");
    fprintf(fp, "Usage:   bcftools [--version|--version-only] [--help] <command> <argument>\n");
    fprintf(fp, "\n");
    fprintf(fp, "Commands:\n");

    int i = 0;
    const char *sep = NULL;
    while (cmds[i].alias)
    {
        if ( !cmds[i].func ) sep = cmds[i].alias;
        if ( sep )
        {
            fprintf(fp, "\n -- %s\n", sep);
            sep = NULL;
        }
        if ( cmds[i].func && cmds[i].help[0]!='-' ) fprintf(fp, "    %-12s %s\n", cmds[i].alias, cmds[i].help);
        i++;
    }
#if ENABLE_BCF_PLUGINS
    fprintf(fp,"\n -- Plugins (collection of programs for calling, file manipulation & analysis)\n");
    int nplugins = count_plugins();
    if ( nplugins )
        fprintf(fp,"    %d plugins available, run \"bcftools plugin -lv\" to see a complete list\n", nplugins);
    else
        fprintf(fp,"    0 plugins available, run \"bcftools plugin -l\" for help\n");
#endif
    fprintf(fp,"\n");
    fprintf(fp,
            " Most commands accept VCF, bgzipped VCF, and BCF with the file type detected\n"
            " automatically even when streaming from a pipe. Indexed VCF and BCF will work\n"
            " in all situations. Un-indexed VCF and BCF and streams will work in most but\n"
            " not all situations.\n");
    fprintf(fp,"\n");
}

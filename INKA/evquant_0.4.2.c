/* This program integrates MaxQuant evidence MSMS Counts
 * for all evidence ids on each line of input file.
 *
 * Written by Alex Henneman
 *
 * Thu Mar 14 10:32:25 CET 2019
 *
 * All rights reserved.
 *
 * Instructions
 * ============
 *
 * Compile as follows:  gcc -o evquant evquant_0.4.2.c
 *
 *  */
#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <getopt.h>
#include  <libgen.h>


/* ------------ */
#define PROGRAM_VERSION_STRING "0.4.2"
#define DEFAULT_EVIDENCE_FILENAME "evidence.txt"
#define DEFAULT_OUTPUT_FILENAME "output.txt"

/* This function should look up in the
 * first line for fields woith the names
 * specified by the constants
 *
 * EV_MSMSCOUNT_TAG
 * EV_EXPERIMENT_TAG */

/* These are the tags we look for in the
 * header of a evidence file */
#define EV_MSMSCOUNT_TAG "MS/MS Count"
#define EV_EXPERIMENT_TAG "Experiment"
#define EV_PSTY_TAG "Phospho (STY)"
/* The size increment step for the
 * arrays conataining experiment index and counts
 * value. */
#define EVIDENCE_ARRAY_LEN_DELTA 100000
#define ERROR_FAILED -1
/* ------------ */

typedef struct {

    unsigned int   len;
    unsigned int   maxlen;
    char         **string;

} string_list_t;

typedef struct {

    unsigned int  len;
    char         *alphabet;

} char_tree_definition_t;

typedef struct CharTreeNode {

    char_tree_definition_t  *def;
    struct CharTreeNode    **branch;
    /*
    unsigned int             nbranches;
    */
    void                    *payload;

} char_tree_node_t;

#define ALPHABET_CAPITALS "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

typedef struct UintVector {

    unsigned int  len;
    unsigned int  maxlen;
    unsigned int *element;

} uint_vector_t;

/* --------- */

#define APPEND_TO_STRING_LIST_DELTA 500

int append_to_string_list_t_dup (

        string_list_t *slist,
        char          *field
){
    int retval=0;

    if (slist->len >= slist->maxlen){
        slist->maxlen += APPEND_TO_STRING_LIST_DELTA;
        if ((slist->string = realloc(slist->string,
                    (slist->maxlen)*sizeof(char *)))==NULL){
            fprintf(stderr,"Failed to expand data in "\
                    "append_to_string_list_t.\n");
            retval = ERROR_FAILED;
        }
    }
    if(retval) return(retval);

    slist->string[slist->len] = strdup(field);
    slist->len++;

    return(retval);
}

/* ------ */

int char_tree_getindex (
        
        char                   *alphabet,
        char                    symbol,
        unsigned int           *idx_p
        
){
    char *found;
    int   retval = 0;

    if ((found=strchr(alphabet,symbol))==NULL){
        fprintf(stderr,"Index of %c not found in "\
                "char_tree_getindex.\n",symbol);
        retval = ERROR_FAILED;
    }
    else {
        *idx_p = found - alphabet;
    }

    return retval;
}

/* ------ */

int char_tree_node_create (
        
        char_tree_definition_t  *def,
        char_tree_node_t       **node_p
        
){
    char_tree_node_t *node;
    int               retval=0;

    if ((node=malloc(sizeof(char_tree_node_t)))==NULL){
        fprintf(stderr,"FAiled allcating char node in char_tree_node_create.\n");
        retval = ERROR_FAILED;
    }
    else {
        node->def       = def;
        /*
        node->nbranches = 0;
        */
        node->branch    = NULL;
        node->payload   = NULL;
        *node_p = node;
    }

    return retval;
}

/* ---------- */
/* This function is responsible for inserting a new sequence
 * into a tree or updating the corresponding node. */

int char_tree_enter (

        char              *path,
        char_tree_node_t  *node,
        int              (*insertfunc) (char_tree_node_t *, void *),
        void              *payload

){
    unsigned int idx;
    int          retval = 0;

    if (strlen(path)){
        if (char_tree_getindex(node->def->alphabet,path[0],&idx)){
            fprintf(stderr,"Failed determining child index in char_tree_enter.\n");
            retval = ERROR_FAILED;
        }
        else {
            if (node->branch==NULL){
                if ((node->branch=calloc(node->def->len,
                                sizeof(char_tree_node_t *)))==NULL){
                    fprintf(stderr,
                            "Failed allocating node branches in char_tree_enter.\n");
                    retval = ERROR_FAILED;
                }
            }
            if (node->branch[idx]==NULL){
                if(char_tree_node_create(node->def,&(node->branch[idx]))){
                    fprintf(stderr,"Failed allocating child node in char_tree_enter.\n");
                    retval = ERROR_FAILED;
                }
            }
            path++;
            if (char_tree_enter(path,node->branch[idx],insertfunc,payload)){
                retval = ERROR_FAILED;
            }
        }
    }
    else {
        if (insertfunc(node,payload)){
            fprintf(stderr,"Failed updating node payload in char_tree_enter.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

int char_tree_destroy_node (

        char_tree_node_t  *node,
        int              (*destroy_func) (void *)

){
    unsigned int i;
    int          retval=0;

    if ((node->branch)==NULL){
        if (destroy_func(node->payload)){
            fprintf(stderr,"Failed destroying payload in char_tree_destroy_node.\n");
            retval = ERROR_FAILED;
        }
    }
    else {
        for(i=0;i<(node->def->len);i++){
            if (node->branch[i]!=NULL){
                if(char_tree_destroy_node(node->branch[i],destroy_func)){
                    retval = ERROR_FAILED;
                }
            }
        }
    }
    free(node);

    return retval;
}

/* ------ */

/* This is probably an entry point function
 * to destroy un used tree. NOTE that this
 * function should only be used once on three
 * root because it destrys the branching definition
 * structure that is is referenced allover the tree.  */

int char_tree_destroy (

        char_tree_node_t  *root,
        int              (*destroy_func) (void *)

){
    char_tree_definition_t *def;
    int                     retval=0;

    if (root==NULL){
        fprintf(stderr,"NULL root poiunter in char_tree_destroy.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((root->def)==NULL){
            fprintf(stderr,"NULL definition pointer in char_tree_destroy.\n");
            retval = ERROR_FAILED;
        }
        else {
            def = root->def;
            if (char_tree_destroy_node(root,destroy_func)){
                fprintf(stderr,"Failed destroying root tree in char_tree_destroy.\n");
                retval = ERROR_FAILED;
            }
            if (def->alphabet) free(def->alphabet);
            free(def);
        }
    }

    return retval;
}

/* ----------- */
#define LINE_BUFFER_LENGTH_DELTA 10000000
#define LINE_BUFFER_LENGTH_INIT  10000000

int fgets_full_line (

        FILE          *fd,
        char         **linebuffer_p,
        unsigned int  *buflen_p
){
    char         *linebuffer;
    unsigned int  buflen;
    int           retval=0;
    off_t         start_offset;

    if ((*linebuffer_p==NULL)||(*buflen_p==0)){
        buflen = LINE_BUFFER_LENGTH_INIT;
        if ((linebuffer=malloc(buflen*sizeof(char)))==NULL){
            fprintf(stderr,"Failed allocating line buffer in "\
                    "fgets_full_line.\n");
            retval = ERROR_FAILED;
        }
        else {
            *linebuffer_p = linebuffer;
            *buflen_p     = buflen;
        }
    }
    else {
        linebuffer = *linebuffer_p;
        buflen     = *buflen_p;
    }
    if(retval)return(retval);

    start_offset = ftello(fd);
    if (fgets(linebuffer,buflen,fd)==NULL){
        /* This is a quiet error because we have
         * reached the end of the file. */
        retval = ERROR_FAILED;
        return(retval);
    }
    while (strlen(linebuffer) >= (buflen-1)){
        buflen += LINE_BUFFER_LENGTH_DELTA;
        if ((linebuffer=realloc(linebuffer,buflen))==NULL){
            fprintf(stderr,"Failed expanding line buffer "\
                    "in fgets_full_line.\n");
            retval = ERROR_FAILED;
            break;
        }
        else {
            *linebuffer_p = linebuffer;
            *buflen_p   = buflen;
        }
        fseeko(fd,start_offset,SEEK_SET);
        if (fgets(linebuffer,buflen,fd)==NULL){
            fprintf(stderr,"Failed getting line in the "\
                    "expansion loop of fgets_full_line.\n");
            retval = ERROR_FAILED;
            break;
        }
    }

    return(retval);
}

/* ----------- */

int allocate_string_list_t (

        unsigned int    len,
        string_list_t **sl_p

){
    string_list_t *sl;
    unsigned int   i;
    int            retval=0;

    if ((sl=malloc(sizeof(string_list_t)))==NULL){
        fprintf(stderr,"Failed to allocate header in allocate_string_list_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((sl->string=malloc(len*sizeof(char *)))==NULL){
            fprintf(stderr,"Failed to allocate data in allocate_string_list_t.\n");
            retval=ERROR_FAILED;
        }
        else {
            sl->len    = 0;
            sl->maxlen = len;
            *sl_p      = sl;
            for(i=0;i<(sl->maxlen);i++){
                sl->string[i] = NULL;
            }
        }
    }
    return(retval);
}

/* ------------- */

int expand_uint_vector_t (

        uint_vector_t *vec,
        unsigned int   delta

){
    int retval=0;

    if (vec==NULL){
        fprintf(stderr,"NULL header pointer in expand_uint_vector_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        vec->maxlen += delta;
        if ((vec->element=realloc(vec->element,
                   (vec->maxlen)*sizeof(unsigned int)))==NULL){
            fprintf(stderr,"Failed to expand uint_vector_t data.\n");
            retval = ERROR_FAILED;
        }
    }

    return retval;
}

/* ---------- */
#define APPEND_TO_UINT_VECTOR_T_DELTA   100

int append_to_uint_vector_t (
        
        unsigned int value,
        uint_vector_t *clist
        
){
    int retval = 0;

    if (clist->len >= clist->maxlen){
        if (expand_uint_vector_t(clist,APPEND_TO_UINT_VECTOR_T_DELTA)){
            fprintf(stderr,"Failed to expand candidate list.\n");
            retval = ERROR_FAILED;
        }
    }
    clist->element[clist->len] = value;
    clist->len++;

    return retval;
}

/* ---------- */

int zero_int_array (

        int          *array,
        unsigned int  len

){
    unsigned int i;
    int retval=0;

    if (array==NULL){
        fprintf(stderr,"NULL array pointer in zero_int_array.\n");
        retval = ERROR_FAILED;
    }
    else {
        for(i=0;i<len;i++){
            array[i] = 0;
        }
    }

    return retval;
}

/* ---------- */

int allocate_uint_vector_t (
        
        unsigned int    len,
        uint_vector_t **vec_pp
        
){
    int            retval = 0;
    uint_vector_t *vec_p;

    if ((vec_p=malloc(sizeof(uint_vector_t)))==NULL){
        fprintf(stderr,"Failed allocating header in "\
                "allocate_uint_vector_t.\n");
        retval = ERROR_FAILED;
    }
    else {
        if ((vec_p->element=malloc(len*sizeof(unsigned int)))==NULL){
            fprintf(stderr,"Failed allocating data in "\
                    "allocate_uint_vector_t.\n");
            retval = ERROR_FAILED;
        }
        else {
            vec_p->len    = 0;
            vec_p->maxlen = len;
            *vec_pp       = vec_p;
        }
    }

    return retval;
}

/* ----------- */

int chomp_string (

        char *line

){
    unsigned int last;
    int          retval=0;

    if (line==NULL){
        fprintf(stderr,"NULL pointer to string in chomp_string.\n");
        retval = ERROR_FAILED;
    }
    else {
        last = strlen(line)-1;
        if (line[last]=='\n') line[last] = '\0';
    }

    return retval;
}

/* ------------ */

int extract_evidence_target_field_indices (
        
        char         *first_line,
        unsigned int *psty_index_p,
        unsigned int *exp_index_p,
        unsigned int *msms_count_index_p
        
){
    char *start,*end,tmp;
    unsigned int field_counter;
    int retval=0,ready;

    /* Note that we go for the experiment field first
     * and assume that they are to be  fond in this
     * order. */
    ready = 0;
    start = first_line;
    field_counter = 1;
    while( (end=strchr(start,'\t'))!=NULL ){
        tmp = *end;
        *end = '\0';
        if (strcmp(start,EV_PSTY_TAG)==0){
            *psty_index_p = field_counter;
            ready = 1;
            break;
        }
        *end = tmp;
        start = end+1;
        field_counter++;
    }
    *end = tmp; start = end+1;
    field_counter++;
    /* Go for the next field. */
    while( (end=strchr(start,'\t'))!=NULL ){
        tmp = *end;
        *end = '\0';
        if (strcmp(start,EV_EXPERIMENT_TAG)==0){
            *exp_index_p = field_counter;
            ready = 1;
            break;
        }
        *end = tmp;
        start = end+1;
        field_counter++;
    }
    /* We skipped this because of the break
     * in the while loop, so we do it here. */
    *end = tmp;
    start = end+1;
    field_counter++;
    /* Check whether everything went ok */
    if (!ready){
        fprintf(stderr,"Failed finding experiment target " \
                "tag in extract_evidence_target_field_indices.\n");
        retval = ERROR_FAILED;
    }
    /* Now we go after the MSMS.Count field */
    ready = 0;
    while( (end=strchr(start,'\t'))!=NULL ){
        tmp = *end;
        *end = '\0';
        if (strcmp(start,EV_MSMSCOUNT_TAG)==0){
            *msms_count_index_p = field_counter;
            ready = 1;
            break;
        }
        *end = tmp;
        start = end+1;
        field_counter++;
    }
    if (!ready){
        fprintf(stderr,"Failed finding MSMS.Count target "\
                "tag in extract_evidence_target_field_indices.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}

/* ------------ */

/* This function extracts given a data line
 * in buffer fields exp_index and msms_count_index
 * into a experiment name string exp_name_p and 
 * count number count_p. 
 * NOTE That we are not making a copy of the
 * experiment name. */

int extract_evidence_line_data (
        
        char          *buffer,
        unsigned int   psty_index,
        unsigned int   exp_index,
        unsigned int   msms_count_index,
        char         **exp_name_p,
        unsigned int  *count_p
        
){
    char *start,*end,*ready_start;
    int   retval=0,times,tot_times,psty_flag;

    /* This is for in case we do not find anything */
    *count_p = 0;
    /* Intialize some values. */
    psty_flag = 0;
    times     = psty_index;
    tot_times = times;
    start     = buffer;
    while (times){
        if ((end=strchr(start,'\t'))==NULL){
            fprintf (stderr,"Failed finding %d-th tab for experiment "\
                    "name in extract_evidence_line_data.\n",tot_times-times+1);
            retval = ERROR_FAILED;
        }
        else {
            *end        = '\0';
            ready_start = start;
            start       = end+1;
            times--;
        }
    }
    *end = '\0';
    if (sscanf(ready_start,"%d",&psty_flag)!=1){
        fprintf(stderr,"Failed reading psty count in "\
                "extract_evidence_line_data.\n");
        fprintf(stderr,"Failed converting psty count:%s.\n",start);
        retval = ERROR_FAILED;
    }
    /* FIXME Add some checking here this is a quick-fix! */
    times       = exp_index - psty_index;
    tot_times   = times;
    while (times){
        if ((end=strchr(start,'\t'))==NULL){
            fprintf (stderr,"Failed finding %d-th tab for experiment "\
                    "name in extract_evidence_line_data.\n",
                    tot_times-times+1);
            retval = ERROR_FAILED;
        }
        else {
            *end        = '\0';
            ready_start = start;
            start       = end+1;
            times--;
        }
    }
    *exp_name_p =  ready_start;
    /* FIXME Add some checking here this is a quick-fix! */
    times       = msms_count_index - exp_index;
    tot_times   = times;
    while (times){
        if ((end=strchr(start,'\t'))==NULL){
            fprintf (stderr,"Failed finding %d-th next tab for "\
                    "count value in extract_evidence_line_data.\n",
                    tot_times-times+1);
            retval = ERROR_FAILED;
        }
        else {
            *end        = '\0';
            ready_start = start;
            start       = end+1;
            times--;
        }
    }

    if (sscanf(ready_start,"%u",count_p)!=1){
        fprintf(stderr,"Failed reading counts in extract_evidence_line_data.\n");
        fprintf(stderr,"Failed converting msms count:%s.\n",ready_start);
        retval = ERROR_FAILED;
    }

    
    return retval;
}

/* ------------ */

typedef struct {

    unsigned int index;

} experiment_node_data_t;

/* ------------ */

typedef struct {

    char          *name;
    string_list_t *list;
    unsigned int   index;

} exp_name_insert_data_t;

/* ------------ */

int insert_exp_name_func (

        char_tree_node_t *node, 
        void             *payload
        
){
    int                     retval=0;
    unsigned int            next_index;
    experiment_node_data_t *data;
    exp_name_insert_data_t *data_pass;

    data_pass = (exp_name_insert_data_t *) payload;

    if (node->payload == NULL){
        if ((data=malloc(sizeof(experiment_node_data_t)))==NULL){
            fprintf(stderr,"Failed allocating new node data in "\
                    "insert_exp_name_func.\n");
            retval = ERROR_FAILED;
        }
        else {
            node->payload = data;
            next_index    = data_pass->list->len;
            /* NOTE That we are dealing with
             * volatile experiment names so we
             * copy them here. */
            if(append_to_string_list_t_dup(data_pass->list,data_pass->name)){
                fprintf(stderr,"Failed appending new experiment "\
                        "name in insert_exp_name_func.\n");
                retval = ERROR_FAILED;
            }
            else {
                /* And we allocate the object
                 * to store the index at the
                 * node */
                if ((data=malloc(sizeof(experiment_node_data_t)))==NULL){
                    fprintf(stderr,"Failed allocatind data container "\
                            "for new end node.\n");
                    retval = ERROR_FAILED;
                }
                else {
                    /* Here we store the C-index of the
                     * name in the list. */
                    data->index   = next_index;
                    node->payload = data;
                    /* And we pass it also backwards
                     * so that it can be used from
                     * the calling function. */
                    data_pass->index = next_index;
                }
            }
        }
    }
    else {
        /* When it's an existing entry
         * we only need to pass back the
         * corresponding stored index. */
        data_pass->index = ( (experiment_node_data_t *) node->payload)->index;
    }

    return retval;
}

/* ------------ */

/* This function extracts the C-array index
 * in the string_list_t of experiment names. */

int translate_experiment_name_to_index (
        
        char             *exp_name,
        char_tree_node_t *root,
        string_list_t    *exp_names,
        unsigned int     *exp_name_index_p
        
){
    exp_name_insert_data_t insert_data;
    int                    retval=0;

    /* Here we fill in the values we
     * are going to need downstream. */
    insert_data.name = exp_name;
    insert_data.list = exp_names;
    /* Make the call into the tree */
    if (char_tree_enter(exp_name,root,insert_exp_name_func,&insert_data)){
        fprintf(stderr,"Failed inserting experiment string "\
                "in translate_experiment_name_to_index.\n");
        retval = ERROR_FAILED;
    }
    else {
        /* When we come out of the tree we have either
         * found our experiment entry before and got
         * its index in insert_data or created a new entry
         * in exp_names and the new index is returned in
         * insert_data.index. */
        *exp_name_index_p = insert_data.index;
    }

    return retval;
}

/* ------------ */

int node_data_destroy_func (
        
        void *data
        
){
    int retval=0;

    if (data==NULL){
        fprintf(stderr,"NULL data pointer in node_data_destroy_func.\n");
        retval = ERROR_FAILED;
    }
    else {
        /* Now this ia struct containing an integer
         * so we just dump the whole struct. */
        free(data);
    }

    return retval;
}

/* ------------ */

int store_evidence_line_data (
        
        unsigned int   exp_index,
        unsigned int   count,
        unsigned int **experiments_p,
        unsigned int **counts_p,
        unsigned int  *evlen_p,
        unsigned int  *ev_maxlen_p
        
){
    unsigned int *e_array,*c_array,len,maxlen;
    int           retval=0;

    len     = *evlen_p;
    maxlen  = *ev_maxlen_p;
    e_array = *experiments_p;
    c_array = *counts_p;
    /* Now we can get to work */
    if (len >= maxlen){
        maxlen += EVIDENCE_ARRAY_LEN_DELTA;
        *ev_maxlen_p = maxlen;
        if ((e_array = realloc(e_array,maxlen*sizeof(int)))==NULL){
            fprintf(stderr,"Failed expanding experiment index "\
                    "array in store_evidence_line_data.\n");
            retval = ERROR_FAILED;
        }
        else {
            *experiments_p = e_array;
        }
        if ((c_array = realloc(c_array,maxlen*sizeof(int)))==NULL){
            fprintf(stderr,"Failed expanding counts array in "\
                    "store_evidence_line_data.\n");
            retval = ERROR_FAILED;
        }
        else {
            *counts_p = c_array;
        }
    }
    /* Here we just add the new data */
    e_array[len] = exp_index;
    c_array[len] = count;
    len++;
    *evlen_p = len;

    return retval;
}

/* ------------ */

int char_tree_init (

        char              *alph,
        char_tree_node_t **root_p

){
    char_tree_node_t       *root;
    char_tree_definition_t *def;
    int                     retval=0;

    if (alph==NULL) {
        fprintf(stderr,"NULL aplhabet pointer in char_tree_init.\n");
        retval = ERROR_FAILED;
    }
    if ((def=malloc(
                sizeof(char_tree_definition_t)))==NULL){
        fprintf(stderr,"Failed allocating alphabet "\
                "structure in char_tree_init.\n");
        retval = ERROR_FAILED;
    }

    if (retval)return retval;

    if (char_tree_node_create(def,&root)){
        fprintf(stderr,"Failed creating root node in char_tree_init.\n");
        retval = ERROR_FAILED;
    }
    else {
        root->def           = def;
        root->def->alphabet =  strdup(alph);
        root->def->len      =  strlen(root->def->alphabet);
        *root_p = root;
    }

    return retval;
}

/* ------------ */
/* This is the initial length
 * of the evidence experiment name
 * index and MSMS count arrays */
#define EVIDENCE_ARRAY_LEN_INIT 100000
/* This is the initial length of the experiment
 * name table length. This table contains
 * all different names. */
#define EXPERIMENT_NAME_LIST_INIT 10

/* These are the characters that are allowed in
 * the experiment names in the evidence table */
#define EXPERIMENT_ALPHABET \
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz_1234567890-+().\\/%"

int load_evidence_file (

        char           *ev_filename,
        char           *alphabet,
        string_list_t **exp_names_p,
        unsigned int  **experiments_p,
        unsigned int  **counts_p,
        unsigned int   *evlen_p

){
    FILE             *ev_fd;
    char             *buffer,*exp_name;
    string_list_t    *experiment_list;
    unsigned int      ev_maxlen,buffer_len,exp_field_index,count,
                      msms_count_field_index,exp_list_index,psty_field_index;
    char_tree_node_t *root;
    int               retval=0;

    if ((ev_fd=fopen(ev_filename,"r"))==NULL){
        fprintf(stderr,"Failed opening evidence file %s.\n",ev_filename);
        retval = ERROR_FAILED;
    }
    if (allocate_string_list_t(EXPERIMENT_NAME_LIST_INIT,exp_names_p)){
        fprintf(stderr,"Failed allocating an experiment name "\
                "list in load_evidence_file.\n");
        retval = ERROR_FAILED;
    }
    if (char_tree_init(alphabet,&root)){
        fprintf(stderr,"failed initializing char tree in load_evidence_file.\n");
        retval = ERROR_FAILED;
    }
    /* We make an initial allocation for the 
     * evidence data arrays. */
    ev_maxlen = EVIDENCE_ARRAY_LEN_INIT;
    if ((*experiments_p=malloc(ev_maxlen*sizeof(int)))==NULL){
        fprintf(stderr,"Failed allocating initial evidence experiment "\
                "name index array in load_evidence_file.\n");
        retval = ERROR_FAILED;
    }
    if ((*counts_p=malloc(ev_maxlen*sizeof(int)))==NULL){
        fprintf(stderr,"Failed allocating initial evidence MSMS "\
                "counts array in load_evidence_file.\n");
        retval = ERROR_FAILED;
    }
    /* FIXME We also need to allocate a char_root_tree or so */
    if (retval) return retval;
    buffer = NULL;

    unsigned int line_count = 0;
    if (fgets_full_line(ev_fd,&buffer,&buffer_len)){
        fprintf(stderr,"Failed reading first line in load_evidence_file.\n");
        retval = ERROR_FAILED;
    }
    else {
        if (extract_evidence_target_field_indices(buffer,&psty_field_index,
                    &exp_field_index,&msms_count_field_index)){
            fprintf(stderr,"Failed extracting target field "\
                    "indices in load_evidence_file.\n");
            retval = ERROR_FAILED;
        }
        else {
            fprintf(stderr,"Fields found in evidence file "\
                    "%s (pSTY,exp,MSMScnt): %u,%u,%u.\n",
                    ev_filename,psty_field_index,exp_field_index,
                    msms_count_field_index);
            *evlen_p = 0;
            /* Now we start processing line
             * by line. */
            while (!fgets_full_line(ev_fd,&buffer,&buffer_len)){
                line_count++;
                if(extract_evidence_line_data(buffer,psty_field_index,
                            exp_field_index,msms_count_field_index,
                            &exp_name,&count)){
                    fprintf(stderr,"Failed extracting data from line "\
                            "%u in load_evidence_file.\n",line_count);
                    retval = ERROR_FAILED;
                    break;
                }
                else {
                    if(translate_experiment_name_to_index(exp_name,root,
                                *exp_names_p,&exp_list_index)){
                        fprintf(stderr,"Failed translating name to index "\
                                "in load_evidence_file.\n");
                        retval = ERROR_FAILED;
                    }
                    else {
                        if (store_evidence_line_data(exp_list_index,count,
                                    experiments_p,counts_p,evlen_p,&ev_maxlen)){
                            fprintf(stderr,"Failed storing line data "\
                                    "in load_evidence_file.\n");
                            retval = ERROR_FAILED;
                        }
                    }
                }
            }
        }
    }

    /* Time to clean up */
    fclose(ev_fd);
    free(buffer);
    /* FIXME We need to clean up
     * the char_tree */
    if (char_tree_destroy(root,node_data_destroy_func)){
        fprintf(stderr,"Failed destroying char tree in load_evidence_file.\n");
        retval = ERROR_FAILED;
    }
    /* Here we truncate the arrays to their exact data length. */
    if ((*experiments_p=realloc(*experiments_p,(*evlen_p)*sizeof(int)))==NULL){
        fprintf(stderr,"Failed truncating experiment index array in "\
                "load_evidence_file.\n");
        retval = ERROR_FAILED;
    }
    if ((*counts_p=realloc(*counts_p,(*evlen_p)*sizeof(int)))==NULL){
        fprintf(stderr,"Failed truncating counts array in load_evidence_file.\n");
        retval = ERROR_FAILED;
    }

    return retval;
}

/* ------------ */

int print_usage (
        
        char *progname

){
    int retval=0;

    printf("\n");
    printf("Usage: %s [options] input.txt\n",progname);
    printf("\n");
    printf("This program reads-in a MaxQuant evidence file and a second\n");
    printf("MaxQuant file that neeeds to be supplied as an argument. The\n");
    printf("program will extract all evidence ids in each line of the\n");
    printf("input file and integrate all MS/MS Counts of the evidence\n");
    printf("ids for each experiment name.\n");
    printf(" The output is a file with a header with the experiment\n");
    printf("names an on each following line the integrated number of\n");
    printf("counts for each experiment. Each line corresponds to\n");
    printf("the same line nr in the input file. If no evidence ids are\n");
    printf("found on a specific line, this results in a row of zero-\n");
    printf("valued integrated counts.\n");
    printf(" Note that the evidence file is assumed to have index 0\n");
    printf("corresponding to the first data line, index 1 with the\n");
    printf("next and so on. This is as for now the case in Max Quant\n");
    printf("files version 1.4 nd 1.5.\n");
    printf("\n");
    printf("Options:\n");
    printf("\n");
    printf(" -h          Print this stuff.\n");
    printf(" -p          Print char tree node alphabet.\n");
    printf(" -e filname  Set evidence file name to flname. Default name is %s.\n",
            DEFAULT_EVIDENCE_FILENAME);
    printf(" -o outnm    Set output file name equal to outnm. Default is %s.\n",
            DEFAULT_OUTPUT_FILENAME);
    printf("\n");
    printf("\n");

    return retval;
}

/* ------------ */

int write_header_to_output (
        
        FILE          *out_fd,
        string_list_t *exp_names
        
){
    unsigned int i;
    int retval=0;

    for(i=0;i<(exp_names->len);i++){
        if (i) fprintf(out_fd,"\t");
        fprintf(out_fd,"\"%s\"",exp_names->string[i]);
    }
    fprintf(out_fd,"\n");

    return retval;
}

/* ------------ */

/* This tag is earched for in the input file
 * which we choose by default to be R output
 * so spaces are replaced by dots */
#define RMOD_EVIDENCE_ID_FIELD_TAB "Evidence.IDs"
#define ORIG_EVIDENCE_ID_FIELD_TAB "Evidence IDs"

int extract_input_file_field_indices (
        
        char         *buffer,
        char         *search_field,
        unsigned int *ev_field_index_p
        
){
    char *start,*end,tmp;
    unsigned int field_counter;
    int retval=0,found;

    found = 0;
    start = buffer;
    field_counter = 1;
    while((end=strchr(start,'\t'))!=NULL){
        tmp = *end;
        *end = '\0';
        if (strcmp(start,search_field) == 0){
            *end = tmp;
            *ev_field_index_p = field_counter;
            found = 1;
            break;
        }
        *end = tmp;
        start = end+1;
        field_counter++;
    }
    if(!found){
        /* Need to make this one quiet
        fprintf(stderr,"Failed finding evidences tab %s in "\
                "extract_input_file_field_indices.\n",
                search_field);
                */
        retval = ERROR_FAILED;
    }

    return retval;
}

/* ------------ */

int expand_scs_uint_list (

        char          *scs_list,
        uint_vector_t *ulist

){
    unsigned int  idx;
    char         *start,*end,tmp;
    int           retval=0;

    start = scs_list;
    /* We emply the list */
    ulist->len = 0;
    /* Here we start parsing the ;-separated list */
    while(end=strchr(start,';')){
        tmp = *end;
        *end = '\0';
        if(sscanf(start,"%u",&idx)!=1){
            fprintf(stderr,"Failed reading and evidence index "\
                    "in extract_row_evidence_ids.\n");
            retval = ERROR_FAILED;
        }
        else {
            *end = tmp;
            if (append_to_uint_vector_t(idx,ulist)){
                fprintf(stderr,"Failed adding evidence index to list "\
                        "in extract_row_evidence_ids.\n");
                retval = ERROR_FAILED;
            }
        }
        start = end+1;
    }
    /* NOTE I am assuming that there is ALWAYS at least a single
     * evidence entry in modificationSpecificPeptides.txt etc. lines. 
     *
     * There should be at least this one:*/
    if(sscanf(start,"%u",&idx)!=1){
        fprintf(stderr,"Failed reading expected evidence index "\
                "in extract_row_evidence_ids.\n");
        retval = ERROR_FAILED;
    }
    else {
        if (append_to_uint_vector_t(idx,ulist)){
            fprintf(stderr,"Failed adding last evidence index to "\
                    "list in extract_row_evidence_ids.\n");
            retval = ERROR_FAILED;
        }
    }


    return retval;
}

/* ------------ */

/* This function extracts from the current line
 * the evidence ids storen in the evidence-ids
 * field and stores these ids in the list ev_ids. */

int extract_row_evidence_ids (
        
        char          *buffer,
        unsigned int   field_nr,
        uint_vector_t *ev_ids
        
){
    char         *start,*end,*tmp_loc,tmp;
    unsigned int  i,target,idx;
    int           retval=0;

    start = buffer;
    target = field_nr - 1;
    for (i=0;i<target;i++){
        if((end=strchr(start,'\t'))==NULL){
            fprintf(stderr,"Failed finding tab in "\
                    "extract_row_evidence_ids.\n");
            retval = ERROR_FAILED;
            break;
        }
        else {
            start = end+1;
        }
    }
    if ((end=strchr(start,'\t'))==NULL){
       fprintf(stderr,"Failed finding end of evidence field in "\
               "extract_row_evidence_ids.\n");
       retval = ERROR_FAILED;
    }
    else {
        tmp     = *end;
        tmp_loc = end;
        *end = '\0';
    }
    if (expand_scs_uint_list (start,ev_ids)){
        fprintf(stderr,"Failed extracting ids from string in "\
                "extract_row_evidence_ids.\n");
        retval = ERROR_FAILED;
    }
    /* Here we restore the closing
     * \t character. */
    *tmp_loc = tmp;

    return retval;
}

/* ------------ */

int write_integrated_row (
        
        FILE         *out_fd,
        unsigned int *row_counts,
        unsigned int  len
        
){
    unsigned int i;
    int retval=0;

    for(i=0;i<len;i++){
        if (i) fprintf(out_fd,"\t");
        fprintf(out_fd,"%u",row_counts[i]);
    }
    fprintf(out_fd,"\n");

    return retval;
}

/* ------------ */

/* This function generates the values
 * of the integrated counts for all
 * experiments corresponding to the
 * current input file line. */

int integrate_row_evidences (
        
        unsigned int  *row_counts,
        unsigned int   row_len,
        uint_vector_t *ev_ids,
        unsigned int  *exp_indices,
        unsigned int  *msms_counts
        
){
    unsigned int i,exp_index,ev_index;
    int          retval=0;

    /* First we clear the row array */
    if (zero_int_array(row_counts,row_len)){
        fprintf(stderr,"Failed clearing row counts.\n");
        retval = ERROR_FAILED;
    }
    /* Next we start integrating */
    for (i=0;i<(ev_ids->len);i++){
        /* These are evidence ID values
         * that start at 1 */
        ev_index = ev_ids->element[i];
        /* Here we assume that experiment indices
         * start at 0 and can be used as C-array
         * indices. */
        exp_index = exp_indices[ev_index];
        row_counts[exp_index] += msms_counts[ev_index];
    }

    return retval;
}

/* ------------ */
/* The initial length of the evidence
 * id list (per line/entry). */
#define EVIDENCE_ID_LIST_LEN_INIT 1000

int main (

        int    argc,
        char **argv

){
    FILE          *out_fd,*in_fd;
    char          *buffer,*ev_filename,c,*out_filename,*input_file,
                  *alphabet;
    string_list_t *exp_names;
    unsigned int  *row_counts,*exp_indices,*msms_counts,psty_field_index,
                   ev_len,buffer_len,ev_field_index,line_count;
    int            retval=0,append_flag;
    uint_vector_t *ev_ids;

    alphabet = strdup(EXPERIMENT_ALPHABET);
    append_flag = 0;
    fprintf(stderr,"\nThis is %s version %s by Alex Henneman.\n",
            basename(argv[0]),PROGRAM_VERSION_STRING);
    /* Here we set the defaults */
    if((out_filename = strdup(DEFAULT_OUTPUT_FILENAME))==NULL){
        fprintf(stderr,"Failed setting default output file name.\n");
        retval = ERROR_FAILED;
    }
    if((ev_filename = strdup(DEFAULT_EVIDENCE_FILENAME))==NULL){
        fprintf(stderr,"Failed setting default evidence file name.\n");
        retval = ERROR_FAILED;
    }
    /* Here we parse the command-line
     * flags. */
    while((c=getopt(argc,argv,"Ao:he:p"))!=-1){
        switch(c){
            case 'p':
                printf("Alphabet:%s\n",alphabet);
                exit(0);
            case 'A':
                append_flag = 1;
                break;
            case 'h':
                print_usage(basename(argv[0]));
                exit(0);
            case 'o':
                free(out_filename);
                if ((out_filename=strdup(optarg))==NULL){
                    fprintf(stderr,"Failed setting output file name from command-line.\n");
                    retval = ERROR_FAILED;
                }
                break;
            case 'e':
                free(ev_filename);
                if ((ev_filename=strdup(optarg))==NULL){
                    fprintf(stderr,"Failed setting evidence file name from command-line.\n");
                    retval = ERROR_FAILED;
                }
                break;
            default:
                print_usage(basename(argv[0]));
                exit(ERROR_FAILED);
        }
    }
    /* Check whether we still have an 
     * input file left on the command-line. */
    if (argc - optind < 1 ){
        print_usage(basename(argv[0]));
        exit(ERROR_FAILED);
    }
    /* We start opening up stuff to get 
     * into it. */
    if((out_fd=fopen(out_filename,"w+"))==NULL){
        fprintf(stderr,"Failed opening output file.\n");
        retval = ERROR_FAILED;
    }
    input_file = argv[optind];
    if((in_fd=fopen(input_file,"r"))==NULL){
        fprintf(stderr,"Failed opening input file %s.\n",input_file);
        retval = ERROR_FAILED;
    }
    /* Here we extract from the evidence file a list
     * of experiment names (in exp_names), and lists
     * of corresponding experiment indices amd msms counts. */
    if (load_evidence_file(ev_filename,alphabet,&exp_names,&exp_indices,
                &msms_counts,&ev_len)){
        fprintf(stderr,"Failed loading evidence file %s.\n",ev_filename);
        retval = ERROR_FAILED;
    }
    else {
        int i;
        printf("Experiment names: ");
        for(i=0;i<(exp_names->len);i++){
            printf("%s ",exp_names->string[i]);
        }
        printf("\n");
    }
    /* Here we allocate a temporary buffer for all
     * the evidence ids, corresponding to each line
     * parsed of the inout file. */
    if (allocate_uint_vector_t(EVIDENCE_ID_LIST_LEN_INIT,&ev_ids)){
            fprintf(stderr,"Failed allocating an evidence id list.\n");
            retval = ERROR_FAILED;
    }

    if(retval) return retval;

    fprintf(stderr,"Loaded %u evidence entries.\n",ev_len);

    /* Here we process the input file row
     * by row. */

    /* First we read the header of the input file
     * to determine which fields we will be extracting on
     * each line. */
    buffer=NULL;
    if (fgets_full_line(in_fd,&buffer,&buffer_len)){
        fprintf(stderr,"Failed reading first line in input file.\n");
        retval = ERROR_FAILED;
    }
    else {
        /* Here we lookup the column number of the evidence
         * id columns using two different column name strings.
         * One is the original MQ name, the second is the R-ifyied
         * string in which all spaces are replaced by dots i.e.. */
        if (extract_input_file_field_indices(buffer,
                    RMOD_EVIDENCE_ID_FIELD_TAB,&ev_field_index)){
            if (extract_input_file_field_indices(buffer,
                        ORIG_EVIDENCE_ID_FIELD_TAB,&ev_field_index)){
                fprintf(stderr,"Failed finding evidence id field index.\n");
                retval = ERROR_FAILED;
            }

        }
        else {
            fprintf(stderr,"Found index in input file %s (Evidence IDs): %u.\n",
                    input_file,ev_field_index);
        }
    }
    /* This object will contain the counts for the matrix
     * row corresponding to the current input file row. */
    if ((row_counts=malloc((exp_names->len)*sizeof(int)))==NULL){
        fprintf(stderr,"Failed allocating row counts array.\n");
        retval = ERROR_FAILED;
    }
    if(retval) return retval;

    /* We prepare the output file 
     * here. */
    if (append_flag){
        chomp_string(buffer);
        fprintf(out_fd,"%s\t",buffer);
    }
    if (write_header_to_output(out_fd,exp_names)){
        fprintf(stderr,"Failed writing output file header.\n");
        retval = ERROR_FAILED;
    }

    /* Here we start processing the input file
     * line by line. We also write each output matrix
     * row one by one to file. */
    line_count = 0;
    while(!fgets_full_line(in_fd,&buffer,&buffer_len)){
        /* Here we extract the list of evidence ids that
         * correspond to this line in the input file. */
        if (extract_row_evidence_ids(buffer,ev_field_index,ev_ids)){
            fprintf(stderr,"Failed extracting evidence ids.\n");
            retval = ERROR_FAILED;
        }
        else {
            /* Here we sum up all the counts corresponding 
             * to this evidence list. */
            if (integrate_row_evidences(row_counts,exp_names->len,
                        ev_ids,exp_indices,msms_counts)){
                fprintf(stderr,"Failed integrating row evidences\n");
                retval = ERROR_FAILED;
            }
            /* And we write out this matrix row. */
            if(append_flag){
                chomp_string(buffer);
                fprintf(out_fd,"%s\t",buffer);
            }
            if (write_integrated_row(out_fd,row_counts,exp_names->len)){
                fprintf(stderr,"Failed writing integrated line.\n");
                retval = ERROR_FAILED;

            }
        }
        line_count++;
    }
    /* And here we are ready processing the
     * input file. */
    fprintf(stderr,"Processed %u rows in input file %s.\n",
            line_count,input_file);
    fprintf(stderr,"Output of %u experiments written to file: %s.\n",
            exp_names->len,out_filename);
    /* We clean up all program crap. */
    fclose(out_fd);
    fclose(in_fd);

    return retval;
}


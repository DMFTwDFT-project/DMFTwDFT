  I`  Ý   k820309    9          19.0        cSR_                                                                                                          
       ../comms.F90 COMMS              COMMS_SETUP COMMS_SETUP_VARS COMMS_END COMMS_BARRIER COMMS_ARRAY_SPLIT ON_ROOT NUM_NODES MY_NODE_ID ROOT_ID gen@COMMS_BCAST gen@COMMS_SEND gen@COMMS_RECV gen@COMMS_REDUCE gen@COMMS_ALLREDUCE gen@COMMS_GATHERV                      @                              
       DP                      @                              
       IO_ERROR                                                        u #COMMS_BCAST_INT    #COMMS_BCAST_LOGICAL    #COMMS_BCAST_REAL 	   #COMMS_BCAST_CMPLX    #COMMS_BCAST_CHAR    #         @     @X                                                 #ARRAY    #SIZE              
D @                                                     
@ @                                          #         @     @X                                                 #ARRAY    #SIZE              
D @                                                     
@ @                                          #         @     @X                             	                    #ARRAY 
   #SIZE              
D @                              
     
                 
@ @                                          #         @     @X                                                 #ARRAY    #SIZE              
D @                                                    
@ @                                          #         @     @X                                                 #ARRAY    #SIZE              
D @                                                  1           
@ @                                                                                                 u #COMMS_SEND_INT    #COMMS_SEND_LOGICAL    #COMMS_SEND_REAL    #COMMS_SEND_CMPLX    #COMMS_SEND_CHAR "   #         @     @X                                                 #ARRAY    #SIZE    #TO              
D @                                                     
@ @                                                    
@ @                                          #         @     @X                                                 #ARRAY    #SIZE    #TO              
D @                                                     
@ @                                                    
@ @                                          #         @     @X                                                 #ARRAY    #SIZE    #TO              
D @                                   
                 
@ @                                                    
@ @                                          #         @     @X                                                 #ARRAY    #SIZE     #TO !             
D @                                                    
@ @                                                     
@ @                               !           #         @     @X                             "                    #ARRAY #   #SIZE $   #TO %             
D @                             #                     1           
@ @                               $                     
@ @                               %                                                                  u #COMMS_RECV_INT &   #COMMS_RECV_LOGICAL *   #COMMS_RECV_REAL .   #COMMS_RECV_CMPLX 2   #COMMS_RECV_CHAR 6   #         @     @X                             &                    #ARRAY '   #SIZE (   #FROM )             
D @                               '                      
@ @                               (                     
@ @                               )           #         @     @X                             *                    #ARRAY +   #SIZE ,   #FROM -             
D @                               +                      
@ @                               ,                     
@ @                               -           #         @     @X                             .                    #ARRAY /   #SIZE 0   #FROM 1             
D @                              /     
                 
@ @                               0                     
@ @                               1           #         @     @X                             2                    #ARRAY 3   #SIZE 4   #FROM 5             
D @                              3                      
@ @                               4                     
@ @                               5           #         @     @X                             6                    #ARRAY 7   #SIZE 8   #FROM 9             
D @                             7                     1           
@ @                               8                     
@ @                               9                                                                  u #COMMS_REDUCE_REAL :   #COMMS_REDUCE_CMPLX >   #COMMS_REDUCE_CMPLX2 B   #COMMS_REDUCE_CMPLX3 G   #         @     @X                             :                    #ARRAY ;   #SIZE <   #OP =             
D @                              ;     
                 
@ @                               <                     
                                =                    1 #         @     @X                             >                    #ARRAY ?   #SIZE @   #OP A             
D @                              ?                      
@ @                               @                     
                                A                    1 #         @     @X                             B                    #ARRAY C   #SIZE1 D   #SIZE2 E   #OP F             
D @                              C                                  &                   &                                                     
                                  D                     
                                  E                     
                                F                    1 #         @     @X                             G                    #ARRAY H   #SIZE1 I   #SIZE2 J   #SIZE3 K   #OP L             
D @                              H                                  &                   &                   &                                                     
                                  I                     
                                  J                     
                                  K                     
                                L                    1                                                        u #COMMS_ALLREDUCE_REAL M   #COMMS_ALLREDUCE_REAL2 Q   #COMMS_ALLREDUCE_REAL3 V   #COMMS_ALLREDUCE_REAL4 \   #COMMS_ALLREDUCE_CMPLX2 c   #COMMS_ALLREDUCE_CMPLX3 h   #COMMS_ALLREDUCE_CMPLX n   #         @     @X                             M                    #ARRAY N   #SIZE O   #OP P             
D @                              N     
                 
@ @                               O                     
                                P                    1 #         @     @X                             Q                    #ARRAY R   #SIZE1 S   #SIZE2 T   #OP U             
D @                              R                   
               &                   &                                                     
                                  S                     
                                  T                     
                                U                    1 #         @     @X                             V                    #ARRAY W   #SIZE1 X   #SIZE2 Y   #SIZE3 Z   #OP [             
D @                              W                   
               &                   &                   &                                                     
                                  X                     
                                  Y                     
                                  Z                     
                                [                    1 #         @     @X                             \                    #ARRAY ]   #SIZE1 ^   #SIZE2 _   #SIZE3 `   #SIZE4 a   #OP b             
D @                              ]                   
               &                   &                   &                   &                                                     
                                  ^                     
                                  _                     
                                  `                     
                                  a                     
                                b                    1 #         @     @X                             c                    #ARRAY d   #SIZE1 e   #SIZE2 f   #OP g             
D @                              d                                  &                   &                   &                                                     
                                  e                     
                                  f                     
                                g                    1 #         @     @X                             h                    #ARRAY i   #SIZE1 j   #SIZE2 k   #SIZE3 l   #OP m             
D @                              i                                  &                   &                   &                                                     
                                  j                     
                                  k                     
                                  l                     
                                m                    1 #         @     @X                             n                    #ARRAY o   #SIZE p   #OP q             
D @                              o                      
@ @                               p                     
                                q                    1                                                     
   u #COMMS_GATHERV_LOGICAL r   #COMMS_GATHERV_REAL_1 y   #COMMS_GATHERV_REAL_2    #COMMS_GATHERV_REAL_3    #COMMS_GATHERV_REAL_2_3    #COMMS_GATHERV_CMPLX_1    #COMMS_GATHERV_CMPLX_2    #COMMS_GATHERV_CMPLX_3    #COMMS_GATHERV_CMPLX_3_4 £   #COMMS_GATHERV_CMPLX_4 ©   #         @     @X                             r                    #ARRAY s   #LOCALCOUNT t   #ROOTGLOBALARRAY u   #COUNTS v   #DISPLS x             
D @                               s                      
@ @                               t                     
D @                               u                     
@ @                               v                     E   p          5 r w       5 r w                              
@ @                               x                     F   p          5 r w       5 r w                     #         @     @X                             y                    #ARRAY z   #LOCALCOUNT {   #ROOTGLOBALARRAY |   #COUNTS }   #DISPLS ~             
D @                              z                   
 !              &                                                     
@ @                               {                     
D @                              |                   
 "              &                                                    
@ @                               }                     #   p          5 r w       5 r w                              
@ @                               ~                     $   p          5 r w       5 r w                     #         @     @X                                                 #ARRAY    #LOCALCOUNT    #ROOTGLOBALARRAY    #COUNTS    #DISPLS              
D @                                                 
 %              &                   &                                                     
@ @                                                    
D @                                                 
 &              &                   &                                                    
@ @                                                    '   p          5 r w       5 r w                              
@ @                                                    (   p          5 r w       5 r w                     #         @     @X                                                 #ARRAY    #LOCALCOUNT    #ROOTGLOBALARRAY    #COUNTS    #DISPLS              
D @                                                 
 )              &                   &                   &                                                     
@ @                                                    
D @                                                 
 *              &                   &                   &                                                    
@ @                                                    +   p          5 r w       5 r w                              
@ @                                                    ,   p          5 r w       5 r w                     #         @     @X                                                 #ARRAY    #LOCALCOUNT    #ROOTGLOBALARRAY    #COUNTS    #DISPLS              
D @                                                 
 -              &                   &                                                     
@ @                                                    
D @                                                 
 .              &                   &                   &                                                    
@ @                                                    /   p          5 r w       5 r w                              
@ @                                                    0   p          5 r w       5 r w                     #         @     @X                                                 #ARRAY    #LOCALCOUNT    #ROOTGLOBALARRAY    #COUNTS    #DISPLS              
D @                                                  1              &                                                     
@ @                                                    
D @                                                  2              &                                                    
@ @                                                    3   p          5 r w       5 r w                              
@ @                                                    4   p          5 r w       5 r w                     #         @     @X                                                 #ARRAY    #LOCALCOUNT    #ROOTGLOBALARRAY    #COUNTS    #DISPLS              
D @                                                  5              &                   &                                                     
@ @                                                    
D @                                                  6              &                   &                                                    
@ @                                                    7   p          5 r w       5 r w                              
@ @                                                    8   p          5 r w       5 r w                     #         @     @X                                                 #ARRAY    #LOCALCOUNT    #ROOTGLOBALARRAY     #COUNTS ¡   #DISPLS ¢             
D @                                                  9              &                   &                   &                                                     
@ @                                                    
D @                                                   :              &                   &                   &                                                    
@ @                               ¡                     ;   p          5 r w       5 r w                              
@ @                               ¢                     <   p          5 r w       5 r w                     #         @     @X                             £                    #ARRAY ¤   #LOCALCOUNT ¥   #ROOTGLOBALARRAY ¦   #COUNTS §   #DISPLS ¨             
D @                              ¤                    =              &                   &                   &                                                     
@ @                               ¥                     
D @                              ¦                    >              &                   &                   &                   &                                                    
@ @                               §                     ?   p          5 r w       5 r w                              
@ @                               ¨                     @   p          5 r w       5 r w                     #         @     @X                             ©                    #ARRAY ª   #LOCALCOUNT «   #ROOTGLOBALARRAY ¬   #COUNTS ­   #DISPLS ®             
D @                              ª                    A              &                   &                   &                   &                                                     
@ @                               «                     
D @                              ¬                    B              &                   &                   &                   &                                                    
@ @                               ­                     C   p          5 r w       5 r w                              
@ @                               ®                     D   p          5 r w       5 r w                                   @                           ¯                          #MPI_UNWEIGHTED °                @                           °                                  @                           ±                          #MPI_WEIGHTS_EMPTY ²                @                           ²                                  @                           ³                          #MPI_BOTTOM ´   #MPI_IN_PLACE µ   #MPI_STATUS_IGNORE ¶                @                           ´                                 @                           µ                                @                           ¶                                p          p            p                                                @                           ·                          #MPI_STATUSES_IGNORE ¸   #MPI_ERRCODES_IGNORE ¹                @                           ¸                                 p          p          p            p          p                                               @                           ¹                                p          p            p                                                @                           º                          #MPI_ARGVS_NULL »   #MPI_ARGV_NULL ¼   -             @                           »                                 p          p          p            p          p                                  -             @                           ¼                                p          p            p                                             @                                ½                      @@                               w                       @@                               ¾                                                         ¿                                                       0#         @                                   À                     #         @                                  Á                     #         @                                   Â                     #         @                                   Ã                     #         @                                   Ä                   #IO!COMMS_ARRAY_SPLIT%MPIFCMB5 Å   #IO!COMMS_ARRAY_SPLIT%MPIFCMB9 Ç   #IO!COMMS_ARRAY_SPLIT%MPIPRIV1 É   #IO!COMMS_ARRAY_SPLIT%MPIPRIV2 Í   #IO!COMMS_ARRAY_SPLIT%MPIPRIVC Ð   #NUMPOINTS Ó   #COUNTS Ô   #DISPLS Õ                 @                        Å                          #COMMS_ARRAY_SPLIT%MPIFCMB5%MPI_UNWEIGHTED Æ                @                         Æ                                  @                        Ç                          #COMMS_ARRAY_SPLIT%MPIFCMB9%MPI_WEIGHTS_EMPTY È                @                         È                                  @                        É                          #COMMS_ARRAY_SPLIT%MPIPRIV1%MPI_BOTTOM Ê   #COMMS_ARRAY_SPLIT%MPIPRIV1%MPI_IN_PLACE Ë   #COMMS_ARRAY_SPLIT%MPIPRIV1%MPI_STATUS_IGNORE Ì                @                         Ê                                 @                         Ë                                @                         Ì                                p          p            p                                                @                        Í                          #COMMS_ARRAY_SPLIT%MPIPRIV2%MPI_STATUSES_IGNORE Î   #COMMS_ARRAY_SPLIT%MPIPRIV2%MPI_ERRCODES_IGNORE Ï                @                         Î                                 p          p          p            p          p                                               @                         Ï                                p          p            p                                                @                        Ð                          #COMMS_ARRAY_SPLIT%MPIPRIVC%MPI_ARGVS_NULL Ñ   #COMMS_ARRAY_SPLIT%MPIPRIVC%MPI_ARGV_NULL Ò   -             @                         Ñ                                 p          p          p            p          p                                  -             @                         Ò                                p          p            p                                            
  @                               Ó                    D                                 Ô                         p           & p          5 r w   n                                       1       5 r w   n                                      1p         p                                            D                                 Õ                         p           & p          5 r w   n                                       1       5 r w   n                                      1p         p                                                fn#fn    »   á   b   uapp(COMMS      C   J  CONSTANTS    ß  I   J  IO     (  ±       gen@COMMS_BCAST     Ù  ]      COMMS_BCAST_INT &   6  @   a   COMMS_BCAST_INT%ARRAY %   v  @   a   COMMS_BCAST_INT%SIZE $   ¶  ]      COMMS_BCAST_LOGICAL *     @   a   COMMS_BCAST_LOGICAL%ARRAY )   S  @   a   COMMS_BCAST_LOGICAL%SIZE !     ]      COMMS_BCAST_REAL '   ð  @   a   COMMS_BCAST_REAL%ARRAY &   0  @   a   COMMS_BCAST_REAL%SIZE "   p  ]      COMMS_BCAST_CMPLX (   Í  @   a   COMMS_BCAST_CMPLX%ARRAY '     @   a   COMMS_BCAST_CMPLX%SIZE !   M  ]      COMMS_BCAST_CHAR '   ª  L   a   COMMS_BCAST_CHAR%ARRAY &   ö  @   a   COMMS_BCAST_CHAR%SIZE    6  ¬       gen@COMMS_SEND    â  e      COMMS_SEND_INT %   G  @   a   COMMS_SEND_INT%ARRAY $     @   a   COMMS_SEND_INT%SIZE "   Ç  @   a   COMMS_SEND_INT%TO #   	  e      COMMS_SEND_LOGICAL )   l	  @   a   COMMS_SEND_LOGICAL%ARRAY (   ¬	  @   a   COMMS_SEND_LOGICAL%SIZE &   ì	  @   a   COMMS_SEND_LOGICAL%TO     ,
  e      COMMS_SEND_REAL &   
  @   a   COMMS_SEND_REAL%ARRAY %   Ñ
  @   a   COMMS_SEND_REAL%SIZE #     @   a   COMMS_SEND_REAL%TO !   Q  e      COMMS_SEND_CMPLX '   ¶  @   a   COMMS_SEND_CMPLX%ARRAY &   ö  @   a   COMMS_SEND_CMPLX%SIZE $   6  @   a   COMMS_SEND_CMPLX%TO     v  e      COMMS_SEND_CHAR &   Û  L   a   COMMS_SEND_CHAR%ARRAY %   '  @   a   COMMS_SEND_CHAR%SIZE #   g  @   a   COMMS_SEND_CHAR%TO    §  ¬       gen@COMMS_RECV    S  g      COMMS_RECV_INT %   º  @   a   COMMS_RECV_INT%ARRAY $   ú  @   a   COMMS_RECV_INT%SIZE $   :  @   a   COMMS_RECV_INT%FROM #   z  g      COMMS_RECV_LOGICAL )   á  @   a   COMMS_RECV_LOGICAL%ARRAY (   !  @   a   COMMS_RECV_LOGICAL%SIZE (   a  @   a   COMMS_RECV_LOGICAL%FROM     ¡  g      COMMS_RECV_REAL &     @   a   COMMS_RECV_REAL%ARRAY %   H  @   a   COMMS_RECV_REAL%SIZE %     @   a   COMMS_RECV_REAL%FROM !   È  g      COMMS_RECV_CMPLX '   /  @   a   COMMS_RECV_CMPLX%ARRAY &   o  @   a   COMMS_RECV_CMPLX%SIZE &   ¯  @   a   COMMS_RECV_CMPLX%FROM     ï  g      COMMS_RECV_CHAR &   V  L   a   COMMS_RECV_CHAR%ARRAY %   ¢  @   a   COMMS_RECV_CHAR%SIZE %   â  @   a   COMMS_RECV_CHAR%FROM !   "  ¡       gen@COMMS_REDUCE "   Ã  e      COMMS_REDUCE_REAL (   (  @   a   COMMS_REDUCE_REAL%ARRAY '   h  @   a   COMMS_REDUCE_REAL%SIZE %   ¨  L   a   COMMS_REDUCE_REAL%OP #   ô  e      COMMS_REDUCE_CMPLX )   Y  @   a   COMMS_REDUCE_CMPLX%ARRAY (     @   a   COMMS_REDUCE_CMPLX%SIZE &   Ù  L   a   COMMS_REDUCE_CMPLX%OP $   %  q      COMMS_REDUCE_CMPLX2 *     ¤   a   COMMS_REDUCE_CMPLX2%ARRAY *   :  @   a   COMMS_REDUCE_CMPLX2%SIZE1 *   z  @   a   COMMS_REDUCE_CMPLX2%SIZE2 '   º  L   a   COMMS_REDUCE_CMPLX2%OP $     |      COMMS_REDUCE_CMPLX3 *     ¼   a   COMMS_REDUCE_CMPLX3%ARRAY *   >  @   a   COMMS_REDUCE_CMPLX3%SIZE1 *   ~  @   a   COMMS_REDUCE_CMPLX3%SIZE2 *   ¾  @   a   COMMS_REDUCE_CMPLX3%SIZE3 '   þ  L   a   COMMS_REDUCE_CMPLX3%OP $   J  þ       gen@COMMS_ALLREDUCE %   H  e      COMMS_ALLREDUCE_REAL +   ­  @   a   COMMS_ALLREDUCE_REAL%ARRAY *   í  @   a   COMMS_ALLREDUCE_REAL%SIZE (   -  L   a   COMMS_ALLREDUCE_REAL%OP &   y  q      COMMS_ALLREDUCE_REAL2 ,   ê  ¤   a   COMMS_ALLREDUCE_REAL2%ARRAY ,     @   a   COMMS_ALLREDUCE_REAL2%SIZE1 ,   Î  @   a   COMMS_ALLREDUCE_REAL2%SIZE2 )     L   a   COMMS_ALLREDUCE_REAL2%OP &   Z  |      COMMS_ALLREDUCE_REAL3 ,   Ö  ¼   a   COMMS_ALLREDUCE_REAL3%ARRAY ,      @   a   COMMS_ALLREDUCE_REAL3%SIZE1 ,   Ò   @   a   COMMS_ALLREDUCE_REAL3%SIZE2 ,   !  @   a   COMMS_ALLREDUCE_REAL3%SIZE3 )   R!  L   a   COMMS_ALLREDUCE_REAL3%OP &   !        COMMS_ALLREDUCE_REAL4 ,   %"  Ô   a   COMMS_ALLREDUCE_REAL4%ARRAY ,   ù"  @   a   COMMS_ALLREDUCE_REAL4%SIZE1 ,   9#  @   a   COMMS_ALLREDUCE_REAL4%SIZE2 ,   y#  @   a   COMMS_ALLREDUCE_REAL4%SIZE3 ,   ¹#  @   a   COMMS_ALLREDUCE_REAL4%SIZE4 )   ù#  L   a   COMMS_ALLREDUCE_REAL4%OP '   E$  q      COMMS_ALLREDUCE_CMPLX2 -   ¶$  ¼   a   COMMS_ALLREDUCE_CMPLX2%ARRAY -   r%  @   a   COMMS_ALLREDUCE_CMPLX2%SIZE1 -   ²%  @   a   COMMS_ALLREDUCE_CMPLX2%SIZE2 *   ò%  L   a   COMMS_ALLREDUCE_CMPLX2%OP '   >&  |      COMMS_ALLREDUCE_CMPLX3 -   º&  ¼   a   COMMS_ALLREDUCE_CMPLX3%ARRAY -   v'  @   a   COMMS_ALLREDUCE_CMPLX3%SIZE1 -   ¶'  @   a   COMMS_ALLREDUCE_CMPLX3%SIZE2 -   ö'  @   a   COMMS_ALLREDUCE_CMPLX3%SIZE3 *   6(  L   a   COMMS_ALLREDUCE_CMPLX3%OP &   (  e      COMMS_ALLREDUCE_CMPLX ,   ç(  @   a   COMMS_ALLREDUCE_CMPLX%ARRAY +   ')  @   a   COMMS_ALLREDUCE_CMPLX%SIZE )   g)  L   a   COMMS_ALLREDUCE_CMPLX%OP "   ³)  N      gen@COMMS_GATHERV &   +        COMMS_GATHERV_LOGICAL ,   +  @   a   COMMS_GATHERV_LOGICAL%ARRAY 1   Ñ+  @   a   COMMS_GATHERV_LOGICAL%LOCALCOUNT 6   ,  @   a   COMMS_GATHERV_LOGICAL%ROOTGLOBALARRAY -   Q,     a   COMMS_GATHERV_LOGICAL%COUNTS -   å,     a   COMMS_GATHERV_LOGICAL%DISPLS %   y-        COMMS_GATHERV_REAL_1 +   	.     a   COMMS_GATHERV_REAL_1%ARRAY 0   .  @   a   COMMS_GATHERV_REAL_1%LOCALCOUNT 5   Õ.     a   COMMS_GATHERV_REAL_1%ROOTGLOBALARRAY ,   a/     a   COMMS_GATHERV_REAL_1%COUNTS ,   õ/     a   COMMS_GATHERV_REAL_1%DISPLS %   0        COMMS_GATHERV_REAL_2 +   1  ¤   a   COMMS_GATHERV_REAL_2%ARRAY 0   ½1  @   a   COMMS_GATHERV_REAL_2%LOCALCOUNT 5   ý1  ¤   a   COMMS_GATHERV_REAL_2%ROOTGLOBALARRAY ,   ¡2     a   COMMS_GATHERV_REAL_2%COUNTS ,   53     a   COMMS_GATHERV_REAL_2%DISPLS %   É3        COMMS_GATHERV_REAL_3 +   Y4  ¼   a   COMMS_GATHERV_REAL_3%ARRAY 0   5  @   a   COMMS_GATHERV_REAL_3%LOCALCOUNT 5   U5  ¼   a   COMMS_GATHERV_REAL_3%ROOTGLOBALARRAY ,   6     a   COMMS_GATHERV_REAL_3%COUNTS ,   ¥6     a   COMMS_GATHERV_REAL_3%DISPLS '   97        COMMS_GATHERV_REAL_2_3 -   É7  ¤   a   COMMS_GATHERV_REAL_2_3%ARRAY 2   m8  @   a   COMMS_GATHERV_REAL_2_3%LOCALCOUNT 7   ­8  ¼   a   COMMS_GATHERV_REAL_2_3%ROOTGLOBALARRAY .   i9     a   COMMS_GATHERV_REAL_2_3%COUNTS .   ý9     a   COMMS_GATHERV_REAL_2_3%DISPLS &   :        COMMS_GATHERV_CMPLX_1 ,   !;     a   COMMS_GATHERV_CMPLX_1%ARRAY 1   ­;  @   a   COMMS_GATHERV_CMPLX_1%LOCALCOUNT 6   í;     a   COMMS_GATHERV_CMPLX_1%ROOTGLOBALARRAY -   y<     a   COMMS_GATHERV_CMPLX_1%COUNTS -   =     a   COMMS_GATHERV_CMPLX_1%DISPLS &   ¡=        COMMS_GATHERV_CMPLX_2 ,   1>  ¤   a   COMMS_GATHERV_CMPLX_2%ARRAY 1   Õ>  @   a   COMMS_GATHERV_CMPLX_2%LOCALCOUNT 6   ?  ¤   a   COMMS_GATHERV_CMPLX_2%ROOTGLOBALARRAY -   ¹?     a   COMMS_GATHERV_CMPLX_2%COUNTS -   M@     a   COMMS_GATHERV_CMPLX_2%DISPLS &   á@        COMMS_GATHERV_CMPLX_3 ,   qA  ¼   a   COMMS_GATHERV_CMPLX_3%ARRAY 1   -B  @   a   COMMS_GATHERV_CMPLX_3%LOCALCOUNT 6   mB  ¼   a   COMMS_GATHERV_CMPLX_3%ROOTGLOBALARRAY -   )C     a   COMMS_GATHERV_CMPLX_3%COUNTS -   ½C     a   COMMS_GATHERV_CMPLX_3%DISPLS (   QD        COMMS_GATHERV_CMPLX_3_4 .   áD  ¼   a   COMMS_GATHERV_CMPLX_3_4%ARRAY 3   E  @   a   COMMS_GATHERV_CMPLX_3_4%LOCALCOUNT 8   ÝE  Ô   a   COMMS_GATHERV_CMPLX_3_4%ROOTGLOBALARRAY /   ±F     a   COMMS_GATHERV_CMPLX_3_4%COUNTS /   EG     a   COMMS_GATHERV_CMPLX_3_4%DISPLS &   ÙG        COMMS_GATHERV_CMPLX_4 ,   iH  Ô   a   COMMS_GATHERV_CMPLX_4%ARRAY 1   =I  @   a   COMMS_GATHERV_CMPLX_4%LOCALCOUNT 6   }I  Ô   a   COMMS_GATHERV_CMPLX_4%ROOTGLOBALARRAY -   QJ     a   COMMS_GATHERV_CMPLX_4%COUNTS -   åJ     a   COMMS_GATHERV_CMPLX_4%DISPLS    yK  d      COMMS!MPIFCMB5    ÝK  H      MPI_UNWEIGHTED    %L  g      COMMS!MPIFCMB9 "   L  H      MPI_WEIGHTS_EMPTY    ÔL        COMMS!MPIPRIV1    ]M  H      MPI_BOTTOM    ¥M  H      MPI_IN_PLACE "   íM  ¤      MPI_STATUS_IGNORE    N        COMMS!MPIPRIV2 $   O  Ä      MPI_STATUSES_IGNORE $   ×O  ¤      MPI_ERRCODES_IGNORE    {P  w      COMMS!MPIPRIVC    òP  Ä      MPI_ARGVS_NULL    ¶Q  ¤      MPI_ARGV_NULL    ZR  @       ON_ROOT    R  @       NUM_NODES    ÚR  @       MY_NODE_ID    S  q       ROOT_ID    S  H       COMMS_SETUP !   ÓS  H       COMMS_SETUP_VARS    T  H       COMMS_END    cT  H       COMMS_BARRIER "   «T        COMMS_ARRAY_SPLIT :   ÉU       IO!COMMS_ARRAY_SPLIT%MPIFCMB5+IO=MPIFCMB5 L   HV  H     COMMS_ARRAY_SPLIT%MPIFCMB5%MPI_UNWEIGHTED+IO=MPI_UNWEIGHTED :   V       IO!COMMS_ARRAY_SPLIT%MPIFCMB9+IO=MPIFCMB9 R   W  H     COMMS_ARRAY_SPLIT%MPIFCMB9%MPI_WEIGHTS_EMPTY+IO=MPI_WEIGHTS_EMPTY :   ZW  Ú     IO!COMMS_ARRAY_SPLIT%MPIPRIV1+IO=MPIPRIV1 D   4X  H     COMMS_ARRAY_SPLIT%MPIPRIV1%MPI_BOTTOM+IO=MPI_BOTTOM H   |X  H     COMMS_ARRAY_SPLIT%MPIPRIV1%MPI_IN_PLACE+IO=MPI_IN_PLACE R   ÄX  ¤     COMMS_ARRAY_SPLIT%MPIPRIV1%MPI_STATUS_IGNORE+IO=MPI_STATUS_IGNORE :   hY  ¸     IO!COMMS_ARRAY_SPLIT%MPIPRIV2+IO=MPIPRIV2 V    Z  Ä     COMMS_ARRAY_SPLIT%MPIPRIV2%MPI_STATUSES_IGNORE+IO=MPI_STATUSES_IGNORE V   äZ  ¤     COMMS_ARRAY_SPLIT%MPIPRIV2%MPI_ERRCODES_IGNORE+IO=MPI_ERRCODES_IGNORE :   [  ­     IO!COMMS_ARRAY_SPLIT%MPIPRIVC+IO=MPIPRIVC L   5\  Ä     COMMS_ARRAY_SPLIT%MPIPRIVC%MPI_ARGVS_NULL+IO=MPI_ARGVS_NULL J   ù\  ¤     COMMS_ARRAY_SPLIT%MPIPRIVC%MPI_ARGV_NULL+IO=MPI_ARGV_NULL ,   ]  @   a   COMMS_ARRAY_SPLIT%NUMPOINTS )   Ý]  6  a   COMMS_ARRAY_SPLIT%COUNTS )   _  6  a   COMMS_ARRAY_SPLIT%DISPLS 
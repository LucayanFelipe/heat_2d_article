#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PERFIL_MEDIO
#if defined(PERFIL_MEDIO)
  #define NX      1000
  #define NY      1000
  #define N_MIN   5
  #define N_MAX   20
  #define WARMUP_RUNS 1
  #define DT      0.75
  #define STEPS   200000
  #define PERFIL_NOME "Medio"

#elif defined(PERFIL_GRANDE)
  #define NX      3000
  #define NY      3000
  #define N_MIN   5
  #define N_MAX   20
  #define WARMUP_RUNS 1
  #define DT      0.75
  #define STEPS   200000
  #define PERFIL_NOME "Grande"

#elif defined(PERFIL_PERSONALIZADO)
  #define NX      100
  #define NY      100
  #define N_MIN   5
  #define N_MAX   20
  #define WARMUP_RUNS 1
  #define DT      0.75
  #define STEPS   25000
  #define PERFIL_NOME "Personalizado"

#else
  #define NX      256
  #define NY      256
  #define N_MIN   30
  #define N_MAX   100
  #define WARMUP_RUNS 3
  #define DT      0.75
  #define STEPS   25000
  #define PERFIL_NOME "Pequeno"
#endif

#define ALPHA   0.1 // -> quão rápido o calor se espalha pelo material.
//Um alpha alto = material conduz bem (cobre). Um alpha baixo = material isola (madeira).
//o valor 0.1 é moderado – representa algo como concreto ou solo úmido.


//DX E DY: o espaçamento entre os pontos da grade nas direções X e Y.
#define DX      1.0
#define DY      1.0

#define EPSILON       0.01    /* erro relativo do IC 95% para convergencia    */
#define OUTLIER_ZSCORE 2.5    /* runs com |z| > 2.5 sao removidos da media   */

#define VALIDATION_STEPS  200
#define VALIDATION_NX     128
#define VALIDATION_NX     128
#define VALIDATION_NY     128
#define VALIDATION_TMIN   1e-6

// debug visual
#define DEBUG_SNAPSHOTS   10
#define DEBUG_LIMIAR      15
#define DEBUG_JANELA_X    5
#define DEBUG_JANELA_Y    5

// exportacao CSV
#define CSV_DIR           "csv_snapshots"
#define CSV_JANELA_RAIO   150 /* raio em celulas da janela central exportada */


#define RX  (ALPHA * DT / (DX * DX))
#define RY  (ALPHA * DT / (DY * DY))



static double grid_a[NX][NY];
static double grid_b[NX][NY];

static double tempos[N_MAX];

static double (*grid)[NY]     = grid_a;
static double (*new_grid)[NY] = grid_b;


void inicializar() {
    memset(grid_a, 0, sizeof(double) * NX * NY);
    memset(grid_b, 0, sizeof(double) * NX * NY);
    grid     = grid_a;
    new_grid = grid_b;
    grid[NX/2][NY/2] = 100.0;
}

void atualizar() {
    double (*g)[NY]  = grid;
    double (*ng)[NY] = new_grid;

    for (int i = 1; i < NX-1; i++) {
        for (int j = 1; j < NY-1; j++) {
            ng[i][j] = g[i][j]
                + RX * (g[i+1][j] - 2.0*g[i][j] + g[i-1][j])
                + RY * (g[i][j+1] - 2.0*g[i][j] + g[i][j-1]);
        }
    }

    for (int j = 0; j < NY; j++) {
        ng[0][j]    = 0.0;
        ng[NX-1][j] = 0.0;
    }
    for (int i = 0; i < NX; i++) {
        ng[i][0]    = 0.0;
        ng[i][NY-1] = 0.0;
    }

    double (*tmp)[NY] = grid;
    grid     = new_grid;
    new_grid = tmp;
}

double energia_total() {
    double soma = 0.0;
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            soma += grid[i][j];
    return soma;
}

double temp_maxima() {
    double mx = 0.0;
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            if (grid[i][j] > mx) mx = grid[i][j];
    return mx;
}


void exportar_csv(int passo) {
    char caminho[256];
    snprintf(caminho, sizeof(caminho), "%s/heat_step_%07d.csv", CSV_DIR, passo);

    FILE *f = fopen(caminho, "w");
    if (!f) {
        fprintf(stderr, "[AVISO] Nao foi possivel abrir %s para escrita.\n", caminho);
        return;
    }

    int cx = NX / 2, cy = NY / 2;
    int ini_i = cx - CSV_JANELA_RAIO; if (ini_i < 0)  ini_i = 0;
    int fim_i = cx + CSV_JANELA_RAIO; if (fim_i > NX) fim_i = NX;
    int ini_j = cy - CSV_JANELA_RAIO; if (ini_j < 0)  ini_j = 0;
    int fim_j = cy + CSV_JANELA_RAIO; if (fim_j > NY) fim_j = NY;

    fprintf(f, "# heat2d_seq - Perfil: %s\n", PERFIL_NOME);
    fprintf(f, "# passo=%d  t_fisico=%.6f  DT=%.6f\n", passo, passo * DT, DT);
    fprintf(f, "# NX=%d  NY=%d  ALPHA=%.6f  DX=%.6f  DY=%.6f\n", NX, NY, ALPHA, DX, DY);
    fprintf(f, "# RX=%.6f  RY=%.6f  STEPS=%d\n", RX, RY, STEPS);
    fprintf(f, "# janela_central: ini_i=%d  fim_i=%d  ini_j=%d  fim_j=%d\n",
            ini_i, fim_i, ini_j, fim_j);
    fprintf(f, "# Tmax=%.6f  Energia=%.6f\n", temp_maxima(), energia_total());
    fprintf(f, "# Colunas: i, j, x, y, T\n");

    for (int i = ini_i; i < fim_i; i++) {
        for (int j = ini_j; j < fim_j; j++) {
            double x = i * DX;
            double y = j * DY;
            fprintf(f, "%d,%d,%.4f,%.4f,%.8f\n", i, j, x, y, grid[i][j]);
        }
    }

    fclose(f);
    printf("  [CSV] Exportado: %s\n", caminho);
}


void debug_snapshot(int passo, int total) {
    int cx = NX / 2;
    int cy = NY / 2;
    int ini_x, fim_x, ini_y, fim_y;

    if (NX <= DEBUG_LIMIAR && NY <= DEBUG_LIMIAR) {
        ini_x = 0; fim_x = NX;
        ini_y = 0; fim_y = NY;
    } else {
        ini_x = cx - DEBUG_JANELA_X; if (ini_x < 0)  ini_x = 0;
        fim_x = cx + DEBUG_JANELA_X; if (fim_x > NX) fim_x = NX;
        ini_y = cy - DEBUG_JANELA_Y; if (ini_y < 0)  ini_y = 0;
        fim_y = cy + DEBUG_JANELA_Y; if (fim_y > NY) fim_y = NY;
    }

    printf("--- Passo %d / %d", passo, total);
    if (NX > DEBUG_LIMIAR || NY > DEBUG_LIMIAR)
        printf("  [janela central x:%d-%d y:%d-%d]", ini_x, fim_x, ini_y, fim_y);
    printf(" ---\n");
    printf("  Tmax=%.4f C   Energia=%.4f\n\n", temp_maxima(), energia_total());

    for (int i = ini_x; i < fim_x; i++) {
        printf("  ");
        for (int j = ini_y; j < fim_y; j++)
            printf("%7.2f ", grid[i][j]);
        printf("\n");
    }
    printf("\n");
}

void debug_visual() {
    #define N_SNAP_LOOP 8

    printf("========================================\n");
    printf(" DEBUG VISUAL - Perfil: %s\n", PERFIL_NOME);
    printf(" %d snapshots em %d passos\n", DEBUG_SNAPSHOTS, STEPS);
    if (NX > DEBUG_LIMIAR || NY > DEBUG_LIMIAR)
        printf(" Grade %dx%d: janela central %dx%d\n",
               NX, NY, DEBUG_JANELA_X*2, DEBUG_JANELA_Y*2);
    printf("========================================\n\n");

    int snap[N_SNAP_LOOP];
    snap[0] = 1;
    for (int k = 0; k < 6; k++)
        snap[1 + k] = 2 + (int)((double)k * (STEPS - 3) / 7.0 + 0.5);
    snap[7] = STEPS - 1;

    for (int k = 1; k < N_SNAP_LOOP; k++)
        if (snap[k] <= snap[k-1])
            snap[k] = snap[k-1] + 1;

    for (int k = 0; k < N_SNAP_LOOP; k++) {
        if (snap[k] < 1)        snap[k] = 1;
        if (snap[k] > STEPS-1)  snap[k] = STEPS - 1;
    }

    inicializar();
    debug_snapshot(0, STEPS);
    exportar_csv(0);

    int next = 0;
    for (int t = 1; t < STEPS; t++) {
        atualizar();
        if (next < N_SNAP_LOOP && t == snap[next]) {
            debug_snapshot(t, STEPS);
            exportar_csv(t);
            next++;
        }
    }

    debug_snapshot(STEPS, STEPS);
    exportar_csv(STEPS);

    #undef N_SNAP_LOOP

    printf("========================================\n");
    printf(" Fim do debug visual (resultado descartado)\n");
    printf("========================================\n\n");
}


void verificar_estabilidade() {
    double rx   = RX;
    double ry   = RY;
    double soma = rx + ry;

    printf("===== VERIFICACAO DE ESTABILIDADE =====\n");
    printf("  Perfil:            %s\n",        PERFIL_NOME);
    printf("  alpha=%.4f  DT=%.4f\n",          ALPHA, DT);
    printf("  DX=%.4f  DY=%.4f\n",             DX, DY);
    printf("  RX = %.6f\n",                    rx);
    printf("  RY = %.6f\n",                    ry);
    printf("  RX + RY = %.6f  (limite 0.5)\n", soma);
    printf("  Intervalo alvo:    [0.10, 0.20]\n");

    if (soma > 0.5) {
        fprintf(stderr, "\n[FATAL] Esquema instavel! RX+RY=%.6f > 0.5\n", soma);
        exit(EXIT_FAILURE);
    }
    if (soma < 0.10 || soma > 0.20)
        printf("  [AVISO] RX+RY=%.4f fora do intervalo alvo\n", soma);

    double alcance = sqrt(4.0 * ALPHA * DT * STEPS);
    double dom_x   = (NX - 1) * DX;
    double dom_y   = (NY - 1) * DY;

    printf("  Status:            ESTAVEL\n");
    printf("  Alcance fisico:    %.2f unidades (%.1f%% em x, %.1f%% em y)\n",
           alcance,
           alcance / dom_x * 100.0,
           alcance / dom_y * 100.0);
    printf("========================================\n\n");
}

void info_problema() {
    double mem_mb   = (double)(NX * NY * sizeof(double)) / (1024.0 * 1024.0);
    double interior = (double)(NX-2) * (NY-2);

    printf("===== INFORMACOES DO PROBLEMA =====\n");
    printf("  Perfil:              %s\n",       PERFIL_NOME);
    printf("  Grade:               %d x %d\n",  NX, NY);
    printf("  Dominio fisico:      %.1f x %.1f\n", (NX-1)*DX, (NY-1)*DY);
    printf("  Memoria por grid:    %.2f MB\n",  mem_mb);
    printf("  Memoria total:       %.2f MB (2 grids)\n", 2.0*mem_mb);
    printf("  Passos de tempo:     %d\n",        STEPS);
    printf("  Celulas interiores:  %lld\n",      (long long)(NX-2)*(NY-2));
    printf("  FLOPs por passo:     %lld\n",      (long long)(interior)*9LL);
    printf("  FLOPs por execucao:  %.3e\n",      interior*9.0*STEPS);
    printf("====================================\n\n");
}

void validar_numerico() {
    printf("===== VALIDACAO NUMERICA =====\n");

    static double vgrid[VALIDATION_NX][VALIDATION_NY];
    static double vnew[VALIDATION_NX][VALIDATION_NY];

    int Nx = VALIDATION_NX;
    int Ny = VALIDATION_NY;
    int cx = Nx / 2;
    int cy = Ny / 2;

    memset(vgrid, 0, sizeof(vgrid));
    vgrid[cx][cy] = 100.0;

    for (int t = 0; t < VALIDATION_STEPS; t++) {
        for (int i = 1; i < Nx-1; i++)
            for (int j = 1; j < Ny-1; j++)
                vnew[i][j] = vgrid[i][j]
                    + RX * (vgrid[i+1][j] - 2.0*vgrid[i][j] + vgrid[i-1][j])
                    + RY * (vgrid[i][j+1] - 2.0*vgrid[i][j] + vgrid[i][j-1]);

        for (int j = 0; j < Ny; j++) { vnew[0][j]    = 0.0; vnew[Nx-1][j] = 0.0; }
        for (int i = 0; i < Nx; i++) { vnew[i][0]    = 0.0; vnew[i][Ny-1] = 0.0; }
        memcpy(vgrid, vnew, sizeof(vgrid));
    }

    double t_fis  = VALIDATION_STEPS * DT;
    double Q      = 100.0 * DX * DY;
    double denom  = 4.0 * M_PI * ALPHA * t_fis;
    double decay  = 4.0 * ALPHA * t_fis;

    double max_abs       = 0.0;
    double max_rel       = 0.0;
    double sum_l2_num    = 0.0;
    double sum_l2_den    = 0.0;
    int    n_validos     = 0;
    int    n_ignorados   = 0;

    for (int i = 1; i < Nx-1; i++) {
        for (int j = 1; j < Ny-1; j++) {
            double x     = (i - cx) * DX;
            double y     = (j - cy) * DY;
            double r2    = x*x + y*y;
            double T_ref = (Q / denom) * exp(-r2 / decay);

            if (T_ref < VALIDATION_TMIN) { n_ignorados++; continue; }

            double T_num = vgrid[i][j];
            double aerr  = fabs(T_num - T_ref);
            double rerr  = aerr / T_ref;

            if (aerr > max_abs) max_abs = aerr;
            if (rerr > max_rel) max_rel = rerr;

            sum_l2_num += aerr * aerr;
            sum_l2_den += T_ref * T_ref;
            n_validos++;
        }
    }

    double erro_l2 = (sum_l2_den > 0.0) ? sqrt(sum_l2_num / sum_l2_den) : 0.0;

    printf("  Grade de validacao:   %d x %d\n",    Nx, Ny);
    printf("  Passos validados:     %d  (t=%.2f)\n", VALIDATION_STEPS, t_fis);
    printf("  Limiar de exclusao:   T_ref < %.0e C\n", VALIDATION_TMIN);
    printf("  Pontos validos:       %d\n",            n_validos);
    printf("  Pontos ignorados:     %d\n",             n_ignorados);
    printf("  Erro absoluto max:    %.6f C\n",        max_abs);
    printf("  Erro relativo max:    %.4f%%\n",        max_rel * 100.0);
    printf("  Erro L2 normalizado:  %.6f (%.4f%%)\n", erro_l2, erro_l2*100.0);

    const char *status;
    if      (erro_l2 < 0.01) status = "APROVADO  (L2 < 1%%)";
    else if (erro_l2 < 0.05) status = "APROVADO  (L2 < 5%%)";
    else if (erro_l2 < 0.10) status = "ACEITAVEL (L2 < 10%%, efeito de borda esperado)";
    else {
        status = "REPROVADO (L2 >= 10%%)";
        fprintf(stderr, "[AVISO] Validacao numerica falhou.\n");
    }

    printf("  Status:               %s\n", status);
    printf("==============================\n\n");
}

void verificar_energia() {
    printf("===== CONSERVACAO DE ENERGIA =====\n");

    inicializar();
    double E0     = energia_total();
    double Eprev  = E0;
    int violacoes = 0;

    printf("  Energia em t=0:    %.6f\n", E0);

    for (int t = 0; t < 1000; t++) {
        atualizar();
        double Ecurr = energia_total();
        if (Ecurr > Eprev + 1e-9) {
            printf("  [ERRO] Energia aumentou no passo %d: %.8f -> %.8f\n",
                   t, Eprev, Ecurr);
            violacoes++;
            if (violacoes > 5) break;
        }
        Eprev = Ecurr;
    }

    printf("  Energia em t=1000: %.6f\n", Eprev);
    printf("  Decaimento total:  %.4f%%\n", (1.0 - Eprev / E0) * 100.0);
    printf("  Status: %s\n", violacoes == 0
           ? "APROVADO (energia monotonicamente decrescente)"
           : "REPROVADO");
    printf("===================================\n\n");
}

double stat_media(double arr[], int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += arr[i];
    return s / n;
}

double stat_desvio(double arr[], int n, double m) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += (arr[i] - m) * (arr[i] - m);
    return (n > 1) ? sqrt(s / (n - 1)) : 0.0;
}

double stat_mediana(double arr[], int n) {
    double tmp[N_MAX];
    memcpy(tmp, arr, sizeof(double) * n);
    for (int i = 1; i < n; i++) {
        double key = tmp[i];
        int j = i - 1;
        while (j >= 0 && tmp[j] > key) { tmp[j+1] = tmp[j]; j--; }
        tmp[j+1] = key;
    }
    return (n % 2 == 0) ? (tmp[n/2-1] + tmp[n/2]) / 2.0 : tmp[n/2];
}

double stat_min(double arr[], int n) {
    double m = arr[0];
    for (int i = 1; i < n; i++) if (arr[i] < m) m = arr[i];
    return m;
}

double stat_max(double arr[], int n) {
    double m = arr[0];
    for (int i = 1; i < n; i++) if (arr[i] > m) m = arr[i];
    return m;
}


int stat_filtrar_outliers(double entrada[], int n, double saida[]) {
    double m = stat_media(entrada, n);
    double s = stat_desvio(entrada, n, m);

    if (s < 1e-12) {
        memcpy(saida, entrada, sizeof(double) * n);
        return n;
    }

    int k = 0;
    for (int i = 0; i < n; i++) {
        double z = fabs(entrada[i] - m) / s;
        if (z <= OUTLIER_ZSCORE)
            saida[k++] = entrada[i];
        else
            printf("  [OUTLIER removido] Run com %.4f s (z=%.2f)\n", entrada[i], z);
    }
    return k;
}


void aquecimento() {
    printf("Aquecimento: %d rodadas de %d passos...\n", WARMUP_RUNS, STEPS);
    fflush(stdout);

    for (int w = 0; w < WARMUP_RUNS; w++) {
        inicializar();
        double t0 = omp_get_wtime();
        for (int t = 0; t < STEPS; t++)
            atualizar();
        double t1 = omp_get_wtime();
        printf("  Warmup %d/%d: %.3f s (descartado)\n", w+1, WARMUP_RUNS, t1-t0);
        fflush(stdout);
    }
    printf("\n");
}

void benchmark() {
    printf("===== BENCHMARK - Perfil %s =====\n", PERFIL_NOME);
    printf("  N_MIN=%d  N_MAX=%d  epsilon=%.1f%%\n",
           N_MIN, N_MAX, EPSILON * 100.0);
    printf("  Filtro outlier: z > %.1f\n\n", OUTLIER_ZSCORE);

    int    n        = 0;
    double erro_rel = 1.0;

    while (n < N_MAX) {

        inicializar();

        double t0 = omp_get_wtime();
        for (int t = 0; t < STEPS; t++)
            atualizar();
        double t1 = omp_get_wtime();

        tempos[n] = t1 - t0;
        n++;

        if (n < N_MIN) {
            printf("  Run %3d | %.4f s  (acumulando)\n", n, tempos[n-1]);
            fflush(stdout);
            continue;
        }

        double filtrados[N_MAX];
        int    nf = stat_filtrar_outliers(tempos, n, filtrados);

        double m  = stat_media(filtrados, nf);
        double s  = stat_desvio(filtrados, nf, m);
        double ci = 1.96 * s / sqrt((double)nf);
        erro_rel  = (m > 0.0) ? ci / m : 1.0;

        printf("  Run %3d | %.4f s | n_fil=%d | media=%.4f | dp=%.4f"
               " | IC95=[%.4f,%.4f] | err=%.2f%%\n",
               n, tempos[n-1], nf, m, s, m-ci, m+ci, erro_rel*100.0);
        fflush(stdout);

        if (erro_rel < EPSILON) {
            printf("\n  Convergencia atingida em %d runs.\n\n", n);
            break;
        }
    }

    if (erro_rel >= EPSILON)
        printf("\n  [AVISO] Nao convergiu em %d runs (err=%.2f%%).\n"
               "  Resultado ainda valido – aumente N_MAX ou verifique\n"
               "  ruido termico/throttling no RPi.\n\n", N_MAX, erro_rel*100.0);

    double filtrados[N_MAX];
    int    nf  = stat_filtrar_outliers(tempos, n, filtrados);

    double m   = stat_media(filtrados, nf);
    double s   = stat_desvio(filtrados, nf, m);
    double ci  = 1.96 * s / sqrt((double)nf);
    double cv  = (m > 0.0) ? s / m : 0.0;
    double med = stat_mediana(filtrados, nf);
    double mn  = stat_min(filtrados, nf);
    double mx  = stat_max(filtrados, nf);

    /* Metricas brutas (sem filtro) para comparacao */
    double m_bruto  = stat_media(tempos, n);
    double mn_bruto = stat_min(tempos, n);
    double mx_bruto = stat_max(tempos, n);

    double interior = (double)(NX-2) * (NY-2);
    double mlups    = (interior * STEPS) / (m * 1.0e6);
    double bw_gbs   = (interior * 6.0 * sizeof(double) * STEPS) / (m * 1.0e9);
    double gflops   = (interior * 9.0 * STEPS) / (m * 1.0e9);

    printf("===== RESULTADOS =====\n");

    printf("\n  -- Tempo (runs filtrados: %d de %d) --\n", nf, n);
    printf("  Media:               %.6f s\n",   m);
    printf("  Mediana:             %.6f s\n",   med);
    printf("  Desvio padrao:       %.6f s\n",   s);
    printf("  Minimo:              %.6f s\n",   mn);
    printf("  Maximo:              %.6f s\n",   mx);
    printf("  IC 95%%:             [%.6f, %.6f] s\n", m-ci, m+ci);
    printf("  Erro relativo IC:    %.4f%%\n",   erro_rel * 100.0);
    printf("  Coef. de variacao:   %.4f%%\n",   cv * 100.0);
    printf("  Outliers removidos:  %d\n",       n - nf);

    printf("\n  -- Tempo bruto (todos os %d runs) --\n", n);
    printf("  Media bruta:         %.6f s\n",   m_bruto);
    printf("  Minimo bruto:        %.6f s\n",   mn_bruto);
    printf("  Maximo bruto:        %.6f s\n",   mx_bruto);

    printf("\n  -- Desempenho --\n");
    printf("  MLUP/s:              %.2f\n",     mlups);
    printf("  GFLOPs/s:            %.4f\n",     gflops);
    printf("  Bandwidth efetiva:   %.2f GB/s\n", bw_gbs);

    printf("\n  -- Escala --\n");
    printf("  Perfil:              %s\n",        PERFIL_NOME);
    printf("  Grade:               %d x %d\n",  NX, NY);
    printf("  Dominio fisico:      %.1f x %.1f\n", (NX-1)*DX, (NY-1)*DY);
    printf("  DT:                  %.4f\n",      DT);
    printf("  RX + RY:             %.4f\n",      RX + RY);
    printf("  Passos:              %d\n",        STEPS);
    printf("  Alcance fisico:      %.2f unidades\n", sqrt(4.0*ALPHA*DT*STEPS));
    printf("  Celulas interiores:  %lld\n",      (long long)(NX-2)*(NY-2));
    printf("  Atualizacoes total:  %.3e\n",      interior * STEPS);

    printf("\n  -- Parametros de benchmark --\n");
    printf("  N_MIN:               %d\n",        N_MIN);
    printf("  N_MAX:               %d\n",        N_MAX);
    printf("  Epsilon:             %.1f%%\n",    EPSILON * 100.0);
    printf("  Warmup runs:         %d\n",        WARMUP_RUNS);
    printf("  Filtro z-score:      %.1f\n",      OUTLIER_ZSCORE);

    printf("\n======================\n");
}

int main() {

    printf("==============================================\n");
    printf("  heat2d_seq - Difusao de Calor 2D\n");
    printf("  Versao: Sequencial  Perfil: %s\n", PERFIL_NOME);
    printf("  Grade: %d x %d   STEPS: %d   DT: %.4f\n",
           NX, NY, STEPS, DT);
    printf("  Benchmark: N_MIN=%d N_MAX=%d epsilon=%.0f%%"
           " warmup=%d z=%.1f\n",
           N_MIN, N_MAX, EPSILON*100.0, WARMUP_RUNS, OUTLIER_ZSCORE);
    printf("==============================================\n\n");

    if (system("mkdir -p " CSV_DIR) != 0)
        fprintf(stderr, "[AVISO] Nao foi possivel criar o diretorio '%s'.\n", CSV_DIR);

    verificar_estabilidade();
    info_problema();
    validar_numerico();
    verificar_energia();
    debug_visual();
    aquecimento();
    benchmark();

    return EXIT_SUCCESS;
}

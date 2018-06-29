/*
 * Programa de Analise Nodal Modificada no Tempo Compacta
 * Por: Matheus Fernandes Moreno (matheus.moreno@poli.ufrj.br)
 *      Paulo Victor Sant'Anna Pedroso de Lima (pv019@poli.ufrj.br)
 *
 * Baseado no programa de demonstracao de analise nodal modificada compacta
 * do professor Antonio Carlos M. de Queiroz.
 *
 * Universidade Federal do Rio de Janeiro - UFRJ
 * Circuitos Eletricos II - 2018.1
 * Professor: Antonio Carlos M. de Queiroz (acmq@coe.ufrj.br)
 *
 * Versao 1.0:
 *   - Diversos bugs consertados. Programa funcional.
 *
 * Versao 0.6:
 *   - Implementado o transistor PNP.
 *   - Implementado um código adicional que diminui a tolerancia do erro e/ou
 *     randomiza alguns valores caso o sistema seja singular.
 *   - Consertado um bug que nao zerava o contador de iteracoes a cada ponto.
 *
 * Versao 0.5:
 *   - Implementadas as iteracoes pelo metodo de Newton-Raphson.
 *   - Implementado o diodo e o transistor NPN.
 *   - Incluido tempo de calculo.
 *
 * Versao 0.3:
 *   - Implementados os elementos reativos (C e L) pelo metodo theta.
 *   - Algumas correcoes de bugs.
 *
 * Versao 0.2:
 *   - Implementada a analise no tempo e calculo do ponto de operacao.
 *   - Implementadas as fontes variantes no tempo SIN e PULSE.
 *
 * Versao 0.1:
 *   - Implementada a analise DC a partir do codigo base.
 *   - Implementado o transformador ideal.
 */


/* Elementos aceitos:
 *
 * Resistor:             R<nome> <no+> <no-> <resistencia>
 * Indutor:              L<nome> <no+> <no-> <indutancia>
 * Capacitor:            C<nome> <no+> <no-> <capacitancia>
 * Fonte de Tensao:      V<nome> <no+> <no-> <parametros>
 * Fonte de Corrente:    I<nome> <no+> <no-> <parametros>
 * Amp. Operacional:     O<nome> <saida+> <saida-> <entrada+> <entrada->
 * Transf. Ideal:        K<nome> <noa> <nob> <noc> <nod> <n>
 * Amp. de Tensao:       E<nome> <noV+> <noV-> <nov+> <nov-> <Av>
 * Amp. de Corrente:     F<nome> <noI+> <noI-> <noi+> <noi-> <Ai>
 * Transcondutor:        G<nome> <noI+> <noI-> <nov+> <nov-> <Gm>
 * Transresistor:        H<nome> <noV+> <noV-> <noi+> <noi-> <Rm>
 * Diodo:                D<nome> <no+> <no-> <Is> <nVt>
 * Transistor bipolar:   Q<nome> <noc> <nob> <noe> <tipo> <a> <ar> <Isbe> <nVtbe> <Isbc> <nVtbc>
 *
 */


#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS   // Impede que o VS reclame de algumas funcoes
#endif

#define versao "1.0 - 06/2018"

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_LINHA             80          // Comprimento maximo de uma linha do netlist
#define MAX_NOME              11          // Comprimento maximo do nome de um elemento
#define MAX_ELEM             200          // Numero maximo de elementos no circuito
#define MAX_NOS              200          // Numero maximo de nos no circuito

#define TOLG_PADRAO        1e-15          // Tolerancia pivotal inicial
#define TOLG_MIN           1e-40          // Tolerancia pivotal minima
#define C_PO                 1e9          // Capacitor no ponto de operacao
#define L_PO                1e-9          // Indutor no ponto de operacao
#define D_VMAX               0.8          // Tensao no qual o diodo se torna linearizado
#define NR_INICIAL           0.1          // Valores iniciais para o Newton-Raphson
#define MAX_ITER             100          // Maximo de iteracoes do Newton-Raphson
#define MAX_RAND              50          // Maximo de randomizacoes pro Newton-Raphson
#define FAIXA_RAND           100          // Faixa de valores na randomizacao
#define TOLG_NR             1e-6          // Tolerancia de erro no NR
#define E_MIN                1.0          // Valor minimo entre calculo de erro absoluto ou relativo

#define THETA_PADRAO        0.55          // Valor padrao pra theta
#define THETA_MIN           1e-6          // Valor minimo pra theta
#define T_PADRAO            1e-3          // Tempo padrao de simulacao
#define PASSO_PADRAO        1e-5          // Passo padrao
#define P_INT_PADRAO           1          // Quantidade padrao de passos internos
#define PI              acos(-1)

// #define DEBUG


typedef struct elemento {    // Elemento do netlist
  char nome[MAX_NOME];
  double valor, j_t0; // Armazenamos a corrente do tempo anterior pra C
  int a, b, c, d, x, y;

  char id[6]; // Pra identificar o tipo de Q e de fonte
  int ciclos;
  double dc, ampl_1, freq, atraso, amort, phi;
  double ampl_2, subida, descida, ligada, periodo; // Valores exclusivos a PULSE

  double Isbe, nVtbe, Isbc, nVtbc; // Para o diodo, usamos somente Isbe e nVtbe
  double alpha, alpha_r;
} elemento;

typedef int tabela[MAX_NOS + 1];


/* Declaracao de variaveis usadas no programa.
 * E mais simples usarmos globais para nao nos preocuparmos com ponteiros. */

elemento netlist[MAX_ELEM];

int
ne,                          // Elementos
nv,                          // Variaveis
neq,                         // Equacoes
nn,                          // Nos
ponto,                       // Ponto de calculo
qtdePontos,                  // Quantidade de pontos
passosInt = P_INT_PADRAO,    // Passos por ponto
customTran = 0,              // Confere se o usuario passou os parametros do tempo
convergiu = 0,               // Afirma se o NR convergiu ou nao
nroIteracoes = 0,            // Conta o numero de iteracoes
nroRandomizacoes = 0,        // Conta o numero de randomizacoes
u, q, r, s;

tabela C, L;    // Vetores para o algoritmo de compactacao

char
nomeArquivo[MAX_NOME + 10], nomeValores[MAX_NOME + 1],
tipo,
na[MAX_NOME], nb[MAX_NOME], nc[MAX_NOME], nd[MAX_NOME],
lista[MAX_NOS + 1][MAX_NOME + 2],
txt[MAX_LINHA + 1],
*param;

FILE *arquivo, *valores;    // Arquivos com o netlist e t0

long double
Yn[MAX_NOS + 1][MAX_NOS + 2],    // Matriz a ser resolvida
tempo = T_PADRAO,                // Tempo de simulacao
passo = PASSO_PADRAO,            // Tamanho do passo
theta = THETA_PADRAO,            // Valor de theta
t0[MAX_NOS + 1],                 // Resultado em t0 para calcular t0 + /\t
en[MAX_NOS + 1],                 // Resultado anterior de NR
inicio, fim,
tolg = TOLG_PADRAO;


/* Resolucao de sistema de equacoes lineares.
 * Metodo de Gauss-Jordan com condensacao pivotal. */
int
resolverSistema(void) {
  int i, j, l, a;
  double t, p;

  for (i = 1; i <= neq; i++) {
    t = 0.0;
    a = i;

    for (l = i; l <= neq; l++) {
      if (fabs(Yn[l][i])>fabs(t)) {
        a = l;
        t = Yn[l][i];
      }
    }

    if (i != a) {
      for (l = 1; l <= neq + 1; l++) {
        p = Yn[i][l];
        Yn[i][l] = Yn[a][l];
        Yn[a][l] = p;
      }
    }

    if (fabs(t) < tolg) {
      return 1;
    }

    for (j = neq + 1; j > 0; j--) {
      Yn[i][j] /= t;
      p = Yn[i][j];
      if (p != 0)
        for (l = 1; l <= neq; l++) {
          if (l != i)
            Yn[l][j] -= Yn[l][i] * p;
        }
    }
  }
  return 0;
}


/* Rotina que zera/inicializa o sistema. */
void
zerarSistema(void) {
  for (u = 0; u <= neq; u++) {
    for (q = 0; q <= neq + 1; q++)
      Yn[u][q] = 0;
  }
}


/* Rotina que conta os nos e atribui numeros a eles. */
int
numeroNo(char *nome) {
  int i, achou;
  i = 0;
  achou = 0;

  while (!achou && i <= nv)
    if (!(achou = !strcmp(nome, lista[i])))
      i++;

  if (!achou) {
    if (nv == MAX_NOS) {
      printf("(!) ERRO: O programa so aceita ate %d nos.\r\n", nv);
      printf("Pressione qualquer tecla para sair...");
      getchar();
      exit(1);
    }
    nv++;
    strcpy(lista[nv], nome);
    return nv; // Novo no
  }

  return i; // No ja conhecido
}


/* Rotina de controle para que o numero de variaveis nao exceda o de nos. */
void
testarNos(void) {
  if (nv > MAX_NOS) {
    printf("(!) ERRO: As variaveis extra excederam o numero de variaveis permitido (%d).\r\n", MAX_NOS);
    printf("Pressione qualquer tecla para sair...");
    getchar();
    exit(1);
  }
}


/* Rotina que le o arquivo com o netlist e cria o vetor de componentes.
 * Tambem e feita a leitura dos parametros para a analise no tempo. */
void
lerNetlist(void) {
  do {
    ne = 0;
    nv = 0;
    strcpy(lista[0], "0");
    printf("Nome do arquivo com o netlist (ex: mna.net): ");
    scanf("%20s", nomeArquivo);
    while (getchar() != '\n') { /* Limpando o buffer de possiveis \n */ }
    arquivo = fopen(nomeArquivo, "r");
    if (arquivo == 0)
      printf("(!) ERRO: Arquivo %s inexistente.\r\n", nomeArquivo);
  } while (arquivo == 0);

  printf("Lendo netlist:\r\n");
  fgets(txt, MAX_LINHA, arquivo);
  printf("Titulo: %s", txt);

  while (fgets(txt, MAX_LINHA, arquivo)) {
    ne++; // Nao usa o netlist[0]
    if (ne > MAX_ELEM) {
      printf("(!) ERRO: O programa so aceita ate %d elementos.\r\n", MAX_ELEM);
      printf("Pressione qualquer tecla para sair...");
      getchar();
      fclose(arquivo);
      exit(1);
    }

    txt[0] = toupper(txt[0]); // O primeiro caractere da linha descreve a linha
    tipo = txt[0];
    sscanf(txt, "%10s", netlist[ne].nome);
    param = txt + strlen(netlist[ne].nome);

    if (tipo == 'R' || tipo == 'L' || tipo == 'C') {
      sscanf(param, "%10s %10s %Lg", na, nb, &netlist[ne].valor);
      printf("%s %s %s %g\r\n", netlist[ne].nome, na, nb, netlist[ne].valor);
      netlist[ne].a = numeroNo(na);
      netlist[ne].b = numeroNo(nb);
    }
    else if (tipo == 'G' || tipo == 'E' || tipo == 'F' || tipo == 'H' || tipo == 'K') {
      sscanf(param, "%10s %10s %10s %10s %Lg", na, nb, nc, nd, &netlist[ne].valor);
      printf("%s %s %s %s %s %g\r\n", netlist[ne].nome, na, nb, nc, nd, netlist[ne].valor);
      netlist[ne].a = numeroNo(na);
      netlist[ne].b = numeroNo(nb);
      netlist[ne].c = numeroNo(nc);
      netlist[ne].d = numeroNo(nd);
    }
    else if (tipo == 'O') {
      sscanf(param, "%10s %10s %10s %10s", na, nb, nc, nd);
      printf("%s %s %s %s %s\r\n", netlist[ne].nome, na, nb, nc, nd);
      netlist[ne].a = numeroNo(na);
      netlist[ne].b = numeroNo(nb);
      netlist[ne].c = numeroNo(nc);
      netlist[ne].d = numeroNo(nd);
    }
    else if (tipo == 'D') {
      sscanf(param, "%10s %10s %Lg %Lg", na, nb, &netlist[ne].Isbe, &netlist[ne].nVtbe);
      printf("%s %s %s %g %g\r\n", netlist[ne].nome, na, nb, netlist[ne].Isbe, netlist[ne].nVtbe);
      netlist[ne].a = numeroNo(na);
      netlist[ne].b = numeroNo(nb);
    }
    else if (tipo == 'Q') {
      sscanf(param, "%10s %10s %10s %10s %Lg %Lg %Lg %Lg %Lg %Lg", nc, nb, na, &netlist[ne].id, &netlist[ne].alpha,
        &netlist[ne].alpha_r, &netlist[ne].Isbe, &netlist[ne].nVtbe, &netlist[ne].Isbc, &netlist[ne].nVtbc);
      printf("%s %s %s %s %s %g %g %g %g %g %g\r\n", netlist[ne].nome, nc, nb, na, netlist[ne].id, netlist[ne].alpha,
        netlist[ne].alpha_r, netlist[ne].Isbe, netlist[ne].nVtbe, netlist[ne].Isbc, netlist[ne].nVtbc);
      netlist[ne].c = numeroNo(nc);
      netlist[ne].b = numeroNo(nb);
      netlist[ne].a = numeroNo(na);
    }
    else if (tipo == 'I' || tipo == 'V') {
      sscanf(param, "%10s %10s %10s", na, nb, &netlist[ne].id);
      param = param + strlen(na) + strlen(nb) + strlen(netlist[ne].id) + 4;
      if (!strcmp(netlist[ne].id, "DC")) {
        sscanf(param, "%Lg", &netlist[ne].valor);
        printf("%s %s %s %s %g\r\n", netlist[ne].nome, na, nb, netlist[ne].id, netlist[ne].valor);
        netlist[ne].a = numeroNo(na);
        netlist[ne].b = numeroNo(nb);
      }
      else if (!strcmp(netlist[ne].id, "SIN")) {
        sscanf(param, "%Lg %Lg %Lg %Lg %Lg %Lg %lu", &netlist[ne].dc, &netlist[ne].ampl_1, &netlist[ne].freq,
          &netlist[ne].atraso, &netlist[ne].amort, &netlist[ne].phi, &netlist[ne].ciclos);
        printf("%s %s %s %s %Lg %Lg %Lg %Lg %Lg %Lg %lu\r\n", netlist[ne].nome, na, nb, netlist[ne].id,
          netlist[ne].dc, netlist[ne].ampl_1, netlist[ne].freq, netlist[ne].atraso, netlist[ne].amort,
          netlist[ne].phi, netlist[ne].ciclos);
        netlist[ne].a = numeroNo(na);
        netlist[ne].b = numeroNo(nb);
      }
      else if (!strcmp(netlist[ne].id, "PULSE")) {
        sscanf(param, "%Lg %Lg %Lg %Lg %Lg %Lg %Lg %lu", &netlist[ne].ampl_1, &netlist[ne].ampl_2,
          &netlist[ne].atraso, &netlist[ne].subida, &netlist[ne].descida, &netlist[ne].ligada,
          &netlist[ne].periodo, &netlist[ne].ciclos);
        printf("%s %s %s %s %Lg %Lg %Lg %Lg %Lg %Lg %Lg %lu\r\n", netlist[ne].nome, na, nb,
          netlist[ne].id, netlist[ne].ampl_1, netlist[ne].ampl_2, netlist[ne].atraso, netlist[ne].subida,
          netlist[ne].descida, netlist[ne].ligada, netlist[ne].periodo, netlist[ne].ciclos);

        netlist[ne].a = numeroNo(na);
        netlist[ne].b = numeroNo(nb);
      }
      else {
        printf("(!) ERRO: Tipo de fonte desconhecido: %s\r\n", netlist[ne].id);
        printf("Pressione qualquer tecla para sair...");
        getchar();
        fclose(arquivo);
        exit(1);
      }
    }
    else if (tipo == '*') { // Comentario comeca com "*"
      printf("Comentario: %s\r", txt);
      ne--;
    }
    else if (tipo == '.') {
      sscanf(param, "%Lg %Lg %*s %Lg %d", &tempo, &passo, &theta, &passosInt);
      if (theta > 1) {
        printf("(!) ERRO: Parametro theta especificado maior que 1.\r\n");
        printf("Pressione qualquer tecla para sair...");
        getchar();
        fclose(arquivo);
        exit(1);
      }
      else if (theta < THETA_MIN)
        theta = THETA_MIN;
      customTran = 1;
      ne--;
    }
    else {
      printf("(!) ERRO: Elemento desconhecido: %s\r\n", txt);
      printf("Pressione qualquer tecla para sair...");
      getchar();
      fclose(arquivo);
      exit(1);
    }
  }
  fclose(arquivo);

  if (!customTran)
    printf("/!\\ Aviso: nao foram passados valores para a analise no tempo. Serao usados os valores padrao.\r\n");
  printf("Tempo de simulacao: %g s\r\n", tempo);
  printf("Tamanho de Passo: %g s\r\n", passo);
  printf("Teta: %g\r\n", theta);
  printf("Passos internos: %d\r\n", passosInt);
  getchar();
}


/* Rotina de simplificacao do sistema com amp. ops. */
void
somar(int *Q, int a, int b) {
  int i, a1, b1;

  if (Q[a] > Q[b]) {
    a1 = Q[b];
    b1 = Q[a];
  }
  else {
    a1 = Q[a];
    b1 = Q[b];
  }

  if (a1 == b1) {
    printf("(!) ERRO: Circuito invalido - Entradas ou saidas de um amp. op. em curto.\r\n");
    printf("Pressione qualquer tecla para sair...\r\n");
    getchar();
    exit(1);
  }

  for (i = 1; i <= MAX_NOS; i++) {
    if (Q[i] == b1)
      Q[i] = a1;
    if (Q[i] > b1)
      Q[i]--;
  }
}


/* Elementos do programa. */
void
operacional(int na, int nb, int nc, int nd) {
#ifdef DEBUG
  printf("Saida: %d %d; entrada %d %d\r\n", na, nb, nc, nd);
#endif
  somar(L, na, nb);
  somar(C, nc, nd);
}

void
transcondutancia(double gm, int n1, int n2, int n3, int n4) {
  Yn[L[n1]][C[n3]] += gm;
  Yn[L[n2]][C[n4]] += gm;
  Yn[L[n1]][C[n4]] -= gm;
  Yn[L[n2]][C[n3]] -= gm;
}

void
condutancia(double g, int a, int b) {
  transcondutancia(g, a, b, a, b);
}

/* Rotina que corrige o tempo de subida e descida de fontes do tipo PULSE. */
void
correcaoPulse(void) {
  for (u = 1; u <= ne; u++) {
    if (!strcmp(netlist[ne].id, "PULSE")) {
      if (netlist[ne].subida < passo)
        netlist[ne].subida = passo;
      if (netlist[ne].descida < passo)
        netlist[ne].descida = passo;
    }
  }
}

/* Funcao que cacula o valor de uma fonte de tensao/corrente no tempo passo*ponto. */
long double
valorFonte(elemento componente) {
  double val, t;

  switch (componente.id[0]) {
  case 'D':
    val = componente.valor;
    break;

  case 'S':
    t = componente.atraso + (componente.ciclos / componente.freq);
    if (passo*ponto <= componente.atraso) // Antes do atraso
      val = componente.dc + componente.ampl_1*(sin(componente.phi * PI / 180));
    else if (passo*ponto >= t) // Fim dos ciclos
      val = componente.dc + componente.ampl_1*exp(-1 * componente.amort * (t - componente.atraso))
      *(sin(2 * PI * componente.freq * (t - componente.atraso) + (componente.phi * PI / 180)));
    else // Caso geral, depois do atraso e antes do fim
      val = componente.dc + componente.ampl_1*exp(-1 * componente.amort * (passo*ponto - componente.atraso))
      *(sin(2 * PI * componente.freq * (passo*ponto - componente.atraso) + (componente.phi * PI / 180)));
    break;

  case 'P':
    if (passo*ponto < componente.ciclos*componente.periodo) {
      t = fmod(passo*ponto, componente.periodo); // Vemos em que momento de um ciclo a fonte esta
      if (t <= componente.atraso)
        val = componente.ampl_1;
      else if (t <= (componente.atraso + componente.subida))
        val = componente.ampl_1 + ((componente.ampl_2 - componente.ampl_1) / componente.subida)
        * (t - componente.atraso);
      else if (t <= (componente.atraso + componente.subida + componente.ligada))
        val = componente.ampl_2;
      else if (t <= (componente.atraso + componente.subida + componente.ligada + componente.descida))
        val = componente.ampl_2 + ((componente.ampl_1 - componente.ampl_2) / componente.descida)
        * (t - (componente.atraso + componente.subida + componente.ligada));
      else
        val = componente.ampl_1;
    }
    else // Se os ciclos terminaram, a fonte fica na amplitude inicial
      val = componente.ampl_1;
    break;
  }

  return val;
}

void
corrente(double i, int a, int b) {
  Yn[L[a]][neq + 1] -= i;
  Yn[L[b]][neq + 1] += i;
}

void
tensao(double v, int a, int b, int x) {
  transcondutancia(1, 0, x, a, b);
  corrente(v, x, 0);
}

void
ganhoTensao(double av, int a, int b, int c, int d, int x) {
  transcondutancia(1, 0, x, a, b);
  transcondutancia(av, x, 0, c, d);
}

void
ganhoCorrente(double ai, int a, int b, int c, int d, int x) {
  transcondutancia(ai, a, b, x, 0);
  transcondutancia(1, c, d, x, 0);
}

void
transresistencia(double rm, int a, int b, int c, int d, int x, int y) {
  transcondutancia(1, 0, y, a, b);
  transcondutancia(rm, y, 0, x, 0);
  transcondutancia(1, c, d, x, 0);
}

void
transformador(double n, int a, int b, int c, int d, int x) {
  Yn[L[a]][C[x]] -= n;
  Yn[L[b]][C[x]] += n;
  Yn[L[c]][C[x]] += 1;
  Yn[L[d]][C[x]] -= 1;
  Yn[L[x]][C[a]] += n;
  Yn[L[x]][C[b]] -= n;
  Yn[L[x]][C[c]] -= 1;
  Yn[L[x]][C[d]] += 1;
}

void
capacitor(double c, int a, int b, double j) {
  condutancia(c / (theta * passo), a, b);
  corrente(c * (t0[a] - t0[b]) / (theta * passo) + ((1 - theta) * j / theta), b, a);
}

void
indutor(double l, int a, int b, int x) {
  if (ponto) {
    Yn[L[a]][C[x]] += 1;
    Yn[L[b]][C[x]] -= 1;
    Yn[L[x]][C[a]] -= 1;
    Yn[L[x]][C[b]] += 1;
    Yn[L[x]][C[x]] += l / (theta * passo);
    Yn[L[x]][neq + 1] += (((1 - theta) / theta) * (t0[a] - t0[b])) + (t0[x] * l / (theta * passo));
  }
  else {
    Yn[L[a]][C[x]] += 1;
    Yn[L[b]][C[x]] -= 1;
    Yn[L[x]][C[a]] -= 1;
    Yn[L[x]][C[b]] += 1;
    Yn[L[x]][C[x]] += L_PO;
  }
}

double
fonteDiodo(double Is, double nVt, int a, int b) {
  double i0, v0;

  if ((nroIteracoes == 0) && (ponto == 0)) v0 = 0.6;
  else v0 = en[a] - en[b];

  if (v0 >= D_VMAX)
    i0 = (Is * (exp(D_VMAX / nVt) - 1)) - ((Is * exp(D_VMAX / nVt) / nVt) * D_VMAX);
  else
    i0 = (Is * (exp(v0 / nVt) - 1))
    - ((Is * exp(v0 / nVt) / nVt) * v0);

  return i0;
}

double
condutanciaDiodo(double Is, double nVt, int a, int b) {
  double g0, v0;

  if ((nroIteracoes == 0) && (ponto == 0)) v0 = 0.6;
  else v0 = en[a] - en[b];

  if (v0 >= D_VMAX)
    g0 = Is * exp(D_VMAX / nVt) / nVt;
  else
    g0 = Is * exp(v0 / nVt) / nVt;

  return g0;
}

void
transistor(double Isbe, double nVtbe, double Isbc, double nVtbc, double al, double al_r,
  int c, int b, int e, char *id) {
  double gc, ge, ic, ie;

  if (!strcmp(id, "NPN")) {
    ge = condutanciaDiodo(Isbe, nVtbe, b, e);
    ie = fonteDiodo(Isbe, nVtbe, b, e);
    gc = condutanciaDiodo(Isbc, nVtbc, b, c);
    ic = fonteDiodo(Isbc, nVtbc, b, c);

    condutancia(ge, b, e);
    corrente(ie, b, e);
    condutancia(gc, b, c);
    corrente(ic, b, c);
    corrente(al * ie, c, b);
    transcondutancia(al * ge, c, b, b, e);
    corrente(al_r * ic, e, b);
    transcondutancia(al_r * gc, e, b, b, c);
  }
  else {
    ge = condutanciaDiodo(Isbe, nVtbe, e, b);
    ie = fonteDiodo(Isbe, nVtbe, e, b);
    gc = condutanciaDiodo(Isbc, nVtbc, c, b);
    ic = fonteDiodo(Isbc, nVtbc, c, b);

    condutancia(ge, e, b);
    corrente(ie, e, b);
    condutancia(gc, c, b);
    corrente(ic, c, b);
    corrente(al * ie, b, c);
    transcondutancia(al * ge, b, c, e, b);
    corrente(al_r * ic, b, e);
    transcondutancia(al_r * gc, b, e, c, b);
  }
}


/* Essa rotina conta os elementos nao aceitos pela analise nodal simples,
 * simplificando com amp. ops. ou nao, dependendo do elemento. */
void
elementosModificada(void) {
  nn = nv;
  neq = nn;

  for (u = 1; u <= ne; u++) {
    tipo = netlist[u].nome[0];
    if (tipo == 'V' || tipo == 'E') {
      nv++;
      strcpy(lista[nv], "j"); // Tem espaco para mais dois caracteres
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
      operacional(netlist[u].a, netlist[u].b, 0, netlist[u].x);
    }
    else if (tipo == 'F') {
      nv++;
      testarNos();
      strcpy(lista[nv], "j"); // Tem espaco para mais dois caracteres
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
      operacional(netlist[u].x, 0, netlist[u].c, netlist[u].d);
    }
    else if (tipo == 'H') {
      nv = nv + 2;
      testarNos();
      strcpy(lista[nv - 1], "jx"); strcat(lista[nv - 1], netlist[u].nome);
      netlist[u].x = nv - 1;
      strcpy(lista[nv], "jy"); strcat(lista[nv], netlist[u].nome);
      netlist[u].y = nv;
      operacional(netlist[u].a, netlist[u].b, 0, netlist[u].y);
      operacional(netlist[u].x, 0, netlist[u].c, netlist[u].d);
    }
    else if (tipo == 'O') {
      operacional(netlist[u].a, netlist[u].b, netlist[u].c, netlist[u].d);
      neq--;
    }
    else if (tipo == 'K') {
      nv++;
      neq++; // Como queremos calcular a corrente, nao usamos ampops
      testarNos();
      strcpy(lista[nv], "j");
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
    }
    else if (tipo == 'L') {
      nv++;
      neq++;
      testarNos();
      strcpy(lista[nv], "j");
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
    }
  }
}


/* Rotina que lista as variaveis e o netlist interno final. */
void
listarTudo(void) {
  printf("Variaveis internas:\r\n");
  for (u = 0; u <= nv; u++)
    printf("%d -> %s (%d)\r\n", u, lista[u], C[u]);
  getchar();

  printf("Netlist interno final\r\n");
  for (u = 1; u <= ne; u++) {
    tipo = netlist[u].nome[0];
    if (tipo == 'R' || tipo == 'I' || tipo == 'V' || tipo == 'C' || tipo == 'L') {
      printf("%s %d %d %g\r\n", netlist[u].nome, netlist[u].a, netlist[u].b, netlist[u].valor);
    }
    else if (tipo == 'G' || tipo == 'E' || tipo == 'F' || tipo == 'H') {
      printf("%s %d %d %d %d %g\r\n", netlist[u].nome, netlist[u].a, netlist[u].b,
        netlist[u].c, netlist[u].d, netlist[u].valor);
    }
    else if (tipo == 'O') {
      printf("%s %d %d %d %d\r\n", netlist[u].nome, netlist[u].a, netlist[u].b,
        netlist[u].c, netlist[u].d);
    }

    if (tipo == 'V' || tipo == 'E' || tipo == 'F' || tipo == 'O' || tipo == 'K' || tipo == 'L')
      printf("Corrente jx: %d\r\n", netlist[u].x);
    else if (tipo == 'H')
      printf("Correntes jx e jy: %d, %d\r\n", netlist[u].x, netlist[u].y);
  }
  getchar();
}


/* Rotina que calcula a corrente no capacitor para ser usada no proximo passo. */
void
memoriaCapacitor(void) {
  double e_a, e_b;

  for (s = 1; s <= ne; s++) {
    if (netlist[s].nome[0] == 'C') {
      if (!ponto) netlist[s].j_t0 = 0.0;
      else {
        if (C[netlist[s].a] == 0) e_a = 0.0;
        else e_a = Yn[C[netlist[s].a]][neq + 1];

        if (C[netlist[s].b] == 0) e_b = 0.0;
        else e_b = Yn[C[netlist[s].b]][neq + 1];

        netlist[s].j_t0 = (netlist[s].valor / (passo * theta)) * ((e_a - e_b)
          - (t0[netlist[s].a] - t0[netlist[s].b])) - (((1 - theta) / theta) * netlist[s].j_t0);
      }
    }
  }
}


/* Rotina que checa a convergencia do Newton-Raphson. */
void
checarConvergencia(void) {
  convergiu = 1;

  for (r = 1; r <= nv; r++) {
    if (C[r]) {
      if (fabs(Yn[C[r]][neq + 1]) > E_MIN) {
        if (fabs((en[r] - Yn[C[r]][neq + 1]) / Yn[C[r]][neq + 1]) > TOLG_NR)
          convergiu = 0;
      }
      else if (fabs(en[r] - Yn[C[r]][neq + 1]) > TOLG_NR)
        convergiu = 0;

#ifdef DEBUG
      printf("Variavel %d, iteracao %d. Convergencia alcancada? %d\r\n", r, nroIteracoes + 1, convergiu);
      printf("Valor anterior: %Lg; Valor atual: %Lg\r\n", en[r], Yn[C[r]][neq + 1]);
      getchar();
#endif

      en[r] = Yn[C[r]][neq + 1];
    }
  }

  nroIteracoes++;
  if ((nroIteracoes == MAX_ITER) && (!convergiu)) {
    for (r = 0; r <= nv; r++) {
      if (C[r])
        en[r] = rand() % FAIXA_RAND;
    }
    nroRandomizacoes++;
    nroIteracoes = 0;
  }
}


/* Rotina que monta as estampas dos elementos. */
void
montarEstampas(void) {
  for (u = 1; u <= ne; u++) {
    switch (netlist[u].nome[0]) {
    case 'R':
      condutancia(1 / netlist[u].valor, netlist[u].a, netlist[u].b);
      break;
    case 'G':
      transcondutancia(netlist[u].valor, netlist[u].a, netlist[u].b,
        netlist[u].c, netlist[u].d);
      break;
    case 'I':
      corrente(valorFonte(netlist[u]), netlist[u].a, netlist[u].b);
      break;
    case 'V':
      tensao(valorFonte(netlist[u]), netlist[u].a, netlist[u].b, netlist[u].x);
      break;
    case 'E':
      ganhoTensao(netlist[u].valor, netlist[u].a, netlist[u].b, netlist[u].c,
        netlist[u].d, netlist[u].x);
      break;
    case 'F':
      ganhoCorrente(netlist[u].valor, netlist[u].a, netlist[u].b, netlist[u].c,
        netlist[u].d, netlist[u].x);
      break;
    case 'H':
      transresistencia(netlist[u].valor, netlist[u].a, netlist[u].b, netlist[u].c,
        netlist[u].d, netlist[u].x, netlist[u].y);
      break;
    case 'K':
      transformador(netlist[u].valor, netlist[u].a, netlist[u].b, netlist[u].c,
        netlist[u].d, netlist[u].x);
      break;
    case 'C':
      if (!ponto) condutancia(1 / C_PO, netlist[u].a, netlist[u].b);
      else capacitor(netlist[u].valor, netlist[u].a, netlist[u].b, netlist[u].j_t0);
      break;
    case 'L':
      indutor(netlist[u].valor, netlist[u].a, netlist[u].b, netlist[u].x);
      break;
    case 'D':
      condutancia(condutanciaDiodo(netlist[u].Isbe, netlist[u].nVtbe, netlist[u].a, netlist[u].b),
        netlist[u].a, netlist[u].b);
      corrente(fonteDiodo(netlist[u].Isbe, netlist[u].nVtbe, netlist[u].a, netlist[u].b),
        netlist[u].a, netlist[u].b);
      break;
    case 'Q':
      transistor(netlist[u].Isbe, netlist[u].nVtbe, netlist[u].Isbc, netlist[u].nVtbc, netlist[u].alpha,
        netlist[u].alpha_r, netlist[u].c, netlist[u].b, netlist[u].a, netlist[u].id);
      break;
    case 'O':
      break;
    }
  }
}


int
main() {
  system("cls");
  printf("Programa de Analise Nodal Modificada Compacta no Tempo (Versao %s)\r\n", versao);
  printf("Por Matheus F. Moreno e Paulo Victor S. Lima\r\n");
  printf("Codigo base por Antonio Carlos M. de Queiroz\r\n");
  srand((unsigned int)time(NULL));

  for (u = 0; u <= MAX_NOS; u++) { // Inicializacao de tabelas
    C[u] = u;
    L[u] = u;
    t0[u] = 0.0;
  }

  lerNetlist(); // Chamada da rotina que le o netlist
  elementosModificada(); // Processamento de elementos da analise modificada
  listarTudo(); // Listagem de variaveis e elementos

  passo = passo / passosInt;
  qtdePontos = (int)round(tempo / passo);
  correcaoPulse();

  strcpy(nomeValores, nomeArquivo); // Cria o nome do arquivo de t0, <nomeArquivo>.tab
  char *pExtensao = strrchr(nomeValores, '.');
  strcpy(pExtensao, ".tab");

  // Inicio da analise
  printf("O circuito tem %d nos, %d variaveis internas, %d equacoes e %d elementos.\r\n", nn, nv, neq, ne);
  getchar();

  valores = fopen(nomeValores, "w");

  // Escreve as variaveis calculadas na primeira linha do .tab
  fprintf(valores, "t");
  for (r = 1; r <= nv; r++) {
    if (C[r] || r <= nn)
      fprintf(valores, " %s", lista[r]);
    if (r == nv)
      fprintf(valores, "\n");
  }

  inicio = clock();

  // Inicializamos o vetor de valores do Newton-Raphson com NR_INICIAL ou 0
  for (r = 0; r <= nv; r++) {
    if (C[r])
      en[r] = NR_INICIAL;
    else
      en[r] = 0.0;
  }

  // Calculo dos pontos
  for (ponto = 0; ponto <= qtdePontos; ponto++) {
    convergiu = 0;
    nroIteracoes = 0;
    nroRandomizacoes = 0;
    tolg = TOLG_PADRAO;

    while ((!convergiu) && (nroRandomizacoes <= MAX_RAND)) {
#ifdef DEBUG
      printf("Iteracao %d, randomizacao %d, ponto %d.\r\n", nroIteracoes + 1, nroRandomizacoes, ponto);
      getchar();
#endif

      zerarSistema();
      montarEstampas();

      /* Esse while diminui a tolerancia caso o sistema seja singular, porque as vezes
      * a singularidade ocorre por conta de erros numericos ou valores ruins no NR. */
      while (resolverSistema()) {
#ifdef DEBUG
        printf("Possivel sistema singular. tolg diminuida para %Lg.\r\n", tolg);
        getchar();
#endif

        if (tolg > TOLG_MIN)
          tolg *= 1e-1;

        if (tolg < TOLG_MIN) {
#ifdef DEBUG
          printf("Sistema singular nao resolvido. Randomizando valores.\r\n");
          printf("Randomizacao %d.\r\n", nroRandomizacoes);
          getchar();
#endif

          for (r = 0; r <= nv; r++) {
            if (C[r])
              en[r] = rand() % FAIXA_RAND;
            else
              en[r] = 0.0;
          }

          nroRandomizacoes++;
          nroIteracoes = 0;
          tolg = TOLG_PADRAO;
          break;
        }
      } /* while resolverSistema */

      checarConvergencia();
    } /* while ((!convergiu) && (nroRandomizacoes <= MAX_RAND)) */

      // Se o sistema nao convergir em um ponto, saimos do programa
    if ((!convergiu) && (nroRandomizacoes > MAX_RAND)) {
      printf("(!) ERRO: O sistema nao convergiu no ponto %d.\r\n", ponto);
      printf("Pressione qualquer tecla para sair...\r\n");
      getchar();
      fclose(valores);
      exit(1);
    }

    // Resultados registrados apos a convergencia do NR
    memoriaCapacitor();
    for (r = 0; r <= nv; r++) {
      if (C[r])
        t0[r] = Yn[C[r]][neq + 1];
      else
        t0[r] = 0.0;
    }

    // Escrevemos os resultados no arquivo
    if (ponto % passosInt == 0) {
      fprintf(valores, "%.7Lg", ponto*passo);
      for (r = 1; r <= nv; r++) {
        if (C[r] || r <= nn)
          fprintf(valores, " %.6Lg", t0[r]);
        if (r == nv)
          fprintf(valores, "\n");
      }
    }
  }

  fim = clock();

  // Informamos a quantidade de pontos calculados e o tempo de analise
  printf("Pronto. %d pontos calculados internamente; %d foram incluidos na tabela.\r\n",
    ponto - 1, (ponto - 1) / passosInt);
  printf("O programa demorou %.4Lg s para simular o circuito.\r\n", (double)(fim - inicio) / CLOCKS_PER_SEC);
  getchar();
  fclose(valores);

  return 0;
}
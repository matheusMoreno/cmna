/*
 * Programa de Analise Nodal Modificada no Tempo Compacta (com interface grafica)
 * Por: Matheus Fernandes Moreno (matheus.moreno@poli.ufrj.br)
 *      Paulo Victor Sant'Anna Pedroso de Lima (pv019@poli.ufrj.br)
 *
 * Baseado no programa de demonstracao de analise nodal modificada compacta
 * do professor Antonio Carlos M. de Queiroz.
 * Implementado com a WinAPI (biblioteca windows.h).
 * Versao 1.0 (22/06/2018).
 *
 * Universidade Federal do Rio de Janeiro - UFRJ
 * Circuitos Eletricos II - 2018.1
 * Professor: Antonio Carlos M. de Queiroz (acmq@coe.ufrj.br)
 *
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS   // Impede que o VS reclame de algumas funcoes
#endif

#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <string>

#include "resource.h"

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

#define ERRO                 500          // Valor de retorno caso haja algum erro nas funcoes

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
  nomeValores[MAX_PATH],
  tipo,
  na[MAX_NOME], nb[MAX_NOME], nc[MAX_NOME], nd[MAX_NOME],
  lista[MAX_NOS + 1][MAX_NOME + 2],
  txt[MAX_LINHA + 1], aux[MAX_LINHA + 1],
  *param;
char nomeArquivo[MAX_PATH] = "circuito.net";

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

/* String para ser printada caso o usuario queira ver como usar o programa. */
std::string ajuda ="\r\nComo utilizar o programa:\r\n"
  "Inclua na primeira linha do .net o título do seu circuito.\r\n\r\n"
  "Os elementos aceitos pelo programa são (que devem ser incluidos um por linha):\r\n"
  " * Resistor:\t\tR<nome> <nó+> <nó-> <resistência>\r\n"
  " * Indutor:\t\t\tL<nome> <nó+> <nó-> <indutância>\r\n"
  " * Capacitor:\t\tC<nome> <nó+> <nó-> <capacitância>\r\n"
  " * Fonte de Tensao:\t\tV<nome> <nó+> <nó-> <parâmetros>\r\n"
  " * Fonte de Corrente:\tI<nome> <nó+> <nó-> <parâmetros>\r\n"
  " * Amp.Operacional:\t\tO<nome> <saída+> <saída-> <entrada+> <entrada->\r\n"
  " * Transf. Ideal:\t\tK<nome> <nó a> <nó b> <nó c> <nó d> <n>\r\n"
  " * Amp. de Tensao:\t\tE<nome> <nó V+> <nó V-> <nó v+> <nó v-> <Av>\r\n"
  " * Amp. de Corrente:\t\tF<nome> <nó I+> <nó I-> <nó i+> <nó i-> <Ai>\r\n"
  " * Transcondutor:\t\tG<nome> <nó I+> <nó I-> <nó v+> <nó v-> <Gm>\r\n"
  " * Transresistor:\t\tH<nome> <nó V+> <nó V-> <nó i+> <nó i-> <Rm>\r\n"
  " * Diodo:\t\t\tD<nome> <nó+> <nó-> <Is> <nVt>\r\n"
  " * Transistor bipolar:\t\tQ<nome> <nó c> <nó b> <nó e> <tipo> <a> <ar> <Isbe> <nVtbe> <Isbc> <nVtbc>\r\n\r\n"
  "No caso das fontes, os parâmetros adicionais serão:\r\n"
  " * DC <valor>\r\n"
  " * SIN <dc> <amplitude> <frequência Hz> <atraso> <amortecimento> <defasagem graus> <ciclos>\r\n"
  " * PULSE <amplitude 1> <amplitude 2> <atraso> <tempo subida> <tempo descida> <tempo ligada> <período> <ciclos>\r\n\r\n"
  "Para informar os parâmetros da simulação, inclua em uma linha do netlist:\r\n"
  ".TRAN <tempo final> <passo> TETA <teta> <passos por ponto>\r\n\r\n"
  "onde 0 < teta =< 1 é o valor a ser usado no método teta de integração.\r\n";

const char g_szClassName[] = "minhaClasseJanela";

std::string str;


/* Rotina que escreve na tela do programa. */
void
my_printf(HWND hPrint, const char *char_str) {
  std::string std_str = char_str;
  int idx = GetWindowTextLength(hPrint);

  SendMessage(hPrint, EM_SETSEL, idx, idx);
  SendMessage(hPrint, EM_REPLACESEL, 0, (LPARAM)std_str.c_str());
}


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
numeroNo(char *nome, HWND hPrint) {
  int i, achou;
  i = 0;
  achou = 0;


  while (!achou && i <= nv)
    if (!(achou = !strcmp(nome, lista[i])))
      i++;

  if (!achou) {
    if (nv == MAX_NOS) {
      sprintf(aux, "\r\n(!) ERRO: O programa só aceita até %d nós.\r\n", nv);
      my_printf(hPrint, aux);
      return ERRO;
    }
    nv++;
    strcpy(lista[nv], nome);
    return nv; // Novo no
  }

  return i; // No ja conhecido
}


/* Rotina de controle para que o numero de variaveis nao exceda o de nos. */
int
testarNos(HWND hPrint) {
  if (nv > MAX_NOS) {
    sprintf(aux, "\r\n(!) ERRO: As variáveis extra excederam o número de variúveis permitido (%d).\r\n", MAX_NOS);
    my_printf(hPrint, aux);
    return ERRO;
  }

  return 0;
}


/* Rotina que le o arquivo com o netlist e cria o vetor de componentes.
 * Tambem e feita a leitura dos parametros para a analise no tempo. */
int
lerNetlist(HWND hPrint) {
  ne = 0;
  nv = 0;
  strcpy(lista[0], "0");
  arquivo = fopen(nomeArquivo, "r");
  if (arquivo == 0) {
    sprintf(aux, "\r\n(!) ERRO: Arquivo %s inexistente.\r\n", nomeArquivo);
    my_printf(hPrint, aux);
    return ERRO;
  }

  sprintf(aux, "\r\nLendo netlist:\r\n");
  my_printf(hPrint, aux);
  fgets(txt, MAX_LINHA, arquivo);
  sprintf(aux, "\r\nTítulo: %s\r\n", txt);
  my_printf(hPrint, aux);

  while (fgets(txt, MAX_LINHA, arquivo)) {
    ne++; // Nao usa o netlist[0]
    if (ne > MAX_ELEM) {
      sprintf(aux, "\r\n(!) ERRO: O programa só aceita até %d elementos.\r\n", MAX_ELEM);
      my_printf(hPrint, aux);
      fclose(arquivo);
      return ERRO;
    }

    txt[0] = toupper(txt[0]); // O primeiro caractere da linha descreve a linha
    tipo = txt[0];
    sscanf(txt, "%10s", netlist[ne].nome);
    param = txt + strlen(netlist[ne].nome);

    if (tipo == 'R' || tipo == 'L' || tipo == 'C') {
      sscanf(param, "%10s %10s %Lg", na, nb, &netlist[ne].valor);
      sprintf(aux, "%s %s %s %g\r\n", netlist[ne].nome, na, nb, netlist[ne].valor);
      my_printf(hPrint, aux);
      netlist[ne].a = numeroNo(na, hPrint);
      netlist[ne].b = numeroNo(nb, hPrint);

      if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO)) {
        fclose(arquivo);
        return ERRO;
      }
    }
    else if (tipo == 'G' || tipo == 'E' || tipo == 'F' || tipo == 'H' || tipo == 'K') {
      sscanf(param, "%10s %10s %10s %10s %Lg", na, nb, nc, nd, &netlist[ne].valor);
      printf("%s %s %s %s %s %g\r\n", netlist[ne].nome, na, nb, nc, nd, netlist[ne].valor);
      netlist[ne].a = numeroNo(na, hPrint);
      netlist[ne].b = numeroNo(nb, hPrint);
      netlist[ne].c = numeroNo(nc, hPrint);
      netlist[ne].d = numeroNo(nd, hPrint);

      if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO) ||
          (netlist[ne].c == ERRO) || (netlist[ne].d == ERRO)) {
        fclose(arquivo);
        return ERRO;
      }
    }
    else if (tipo == 'O') {
      sscanf(param, "%10s %10s %10s %10s", na, nb, nc, nd);
      printf("%s %s %s %s %s\r\n", netlist[ne].nome, na, nb, nc, nd);
      netlist[ne].a = numeroNo(na, hPrint);
      netlist[ne].b = numeroNo(nb, hPrint);
      netlist[ne].c = numeroNo(nc, hPrint);
      netlist[ne].d = numeroNo(nd, hPrint);

      if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO) ||
        (netlist[ne].c == ERRO) || (netlist[ne].d == ERRO)) {
        fclose(arquivo);
        return ERRO;
      }
    }
    else if (tipo == 'D') {
      sscanf(param, "%10s %10s %Lg %Lg", na, nb, &netlist[ne].Isbe, &netlist[ne].nVtbe);
      printf("%s %s %s %g %g\r\n", netlist[ne].nome, na, nb, netlist[ne].Isbe, netlist[ne].nVtbe);
      netlist[ne].a = numeroNo(na, hPrint);
      netlist[ne].b = numeroNo(nb, hPrint);

      if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO)) {
        fclose(arquivo);
        return ERRO;
      }
    }
    else if (tipo == 'Q') {
      sscanf(param, "%10s %10s %10s %10s %Lg %Lg %Lg %Lg %Lg %Lg", nc, nb, na, &netlist[ne].id, &netlist[ne].alpha,
        &netlist[ne].alpha_r, &netlist[ne].Isbe, &netlist[ne].nVtbe, &netlist[ne].Isbc, &netlist[ne].nVtbc);
      printf("%s %s %s %s %s %g %g %g %g %g %g\r\n", netlist[ne].nome, nc, nb, na, netlist[ne].id, netlist[ne].alpha,
        netlist[ne].alpha_r, netlist[ne].Isbe, netlist[ne].nVtbe, netlist[ne].Isbc, netlist[ne].nVtbc);
      netlist[ne].c = numeroNo(nc, hPrint);
      netlist[ne].b = numeroNo(nb, hPrint);
      netlist[ne].a = numeroNo(na, hPrint);

      if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO) ||
        (netlist[ne].c == ERRO)) {
        fclose(arquivo);
        return ERRO;
      }
    }
    else if (tipo == 'I' || tipo == 'V') {
      sscanf(param, "%10s %10s %10s", na, nb, &netlist[ne].id);
      param = param + strlen(na) + strlen(nb) + strlen(netlist[ne].id) + 4;
      if (!strcmp(netlist[ne].id, "DC")) {
        sscanf(param, "%Lg", &netlist[ne].valor);
        printf("%s %s %s %s %g\r\n", netlist[ne].nome, na, nb, netlist[ne].id, netlist[ne].valor);
        netlist[ne].a = numeroNo(na, hPrint);
        netlist[ne].b = numeroNo(nb, hPrint);

        if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO)) {
          fclose(arquivo);
          return ERRO;
        }
      }
      else if (!strcmp(netlist[ne].id, "SIN")) {
        sscanf(param, "%Lg %Lg %Lg %Lg %Lg %Lg %lu", &netlist[ne].dc, &netlist[ne].ampl_1, &netlist[ne].freq,
          &netlist[ne].atraso, &netlist[ne].amort, &netlist[ne].phi, &netlist[ne].ciclos);
        printf("%s %s %s %s %Lg %Lg %Lg %Lg %Lg %Lg %lu\r\n", netlist[ne].nome, na, nb, netlist[ne].id,
          netlist[ne].dc, netlist[ne].ampl_1, netlist[ne].freq, netlist[ne].atraso, netlist[ne].amort,
          netlist[ne].phi, netlist[ne].ciclos);
        netlist[ne].a = numeroNo(na, hPrint);
        netlist[ne].b = numeroNo(nb, hPrint);

        if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO)) {
          fclose(arquivo);
          return ERRO;
        }
      }
      else if (!strcmp(netlist[ne].id, "PULSE")) {
        sscanf(param, "%Lg %Lg %Lg %Lg %Lg %Lg %Lg %lu", &netlist[ne].ampl_1, &netlist[ne].ampl_2,
          &netlist[ne].atraso, &netlist[ne].subida, &netlist[ne].descida, &netlist[ne].ligada,
          &netlist[ne].periodo, &netlist[ne].ciclos);
        printf("%s %s %s %s %Lg %Lg %Lg %Lg %Lg %Lg %Lg %lu\r\n", netlist[ne].nome, na, nb,
          netlist[ne].id, netlist[ne].ampl_1, netlist[ne].ampl_2, netlist[ne].atraso, netlist[ne].subida,
          netlist[ne].descida, netlist[ne].ligada, netlist[ne].periodo, netlist[ne].ciclos);

        netlist[ne].a = numeroNo(na, hPrint);
        netlist[ne].b = numeroNo(nb, hPrint);

        if ((netlist[ne].a == ERRO) || (netlist[ne].b == ERRO)) {
          fclose(arquivo);
          return ERRO;
        }
      }
      else {
          sprintf(aux, "\r\n(!) ERRO: Tipo de fonte desconhecido: %s\r\n", netlist[ne].id);
          my_printf(hPrint, aux);
          fclose(arquivo);
          return ERRO;
        }
    }
    else if (tipo == '*') { // Comentario comeca com "*"
      sprintf(aux, "Comentário: %s\r", txt);
      my_printf(hPrint, aux);
      fclose(arquivo);
      return ERRO;
      ne--;
    }
    else if (tipo == '.') {
      sscanf(param, "%Lg %Lg %*s %Lg %d", &tempo, &passo, &theta, &passosInt);
      if (theta > 1) {
        my_printf(hPrint, "\r\n(!) ERRO: Parâmetro teta especificado maior que 1.\r\n");
        fclose(arquivo);
        return ERRO;
      }
      else if (theta < THETA_MIN)
        theta = THETA_MIN;
      customTran = 1;
      ne--;
    }
    else {
      sprintf(aux, "\r\n(!) ERRO: Elemento desconhecido: %s\r\n", txt);
      my_printf(hPrint, aux);
      fclose(arquivo);
      return ERRO;
    }
  }
  fclose(arquivo);

  if (!customTran)
    my_printf(hPrint, "/!\\ Aviso: não foram passados valores para a análise no tempo. Serão usados os valores padrão.\r\n");
  sprintf(aux, "Tempo de simulação: %g s\r\n", tempo);
  my_printf(hPrint, aux);
  sprintf(aux, "Tamanho de Passo: %g s\r\n", passo);
  my_printf(hPrint, aux);
  sprintf(aux, "Teta: %g\r\n", theta);
  my_printf(hPrint, aux);
  sprintf(aux, "Passos internos: %d\r\n", passosInt);
  my_printf(hPrint, aux);

  return 0;
}


/* Rotina de simplificacao do sistema com amp. ops. */
int
somar(int *Q, int a, int b, HWND hPrint) {
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
    my_printf(hPrint, "(!) ERRO: Circuito inválido - Entradas ou saídas de um amp. op. em curto.\r\n");
    return ERRO;
  }

  for (i = 1; i <= MAX_NOS; i++) {
    if (Q[i] == b1)
      Q[i] = a1;
    if (Q[i] > b1)
      Q[i]--;
  }

  return 0;
}


/* Elementos do programa. */
int
operacional(int na, int nb, int nc, int nd, HWND hPrint) {
  if (somar(L, na, nb, hPrint)) return ERRO;
  if (somar(C, nc, nd, hPrint)) return ERRO;

  return 0;
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
int
elementosModificada(HWND hPrint) {
  nn = nv;
  neq = nn;

  for (u = 1; u <= ne; u++) {
    tipo = netlist[u].nome[0];
    if (tipo == 'V' || tipo == 'E') {
      nv++;
      strcpy(lista[nv], "j"); // Tem espaco para mais dois caracteres
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
      if (operacional(netlist[u].a, netlist[u].b, 0, netlist[u].x, hPrint)) return ERRO;
    }
    else if (tipo == 'F') {
      nv++;
      if (testarNos(hPrint)) return ERRO;
      strcpy(lista[nv], "j"); // Tem espaco para mais dois caracteres
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
      if (operacional(netlist[u].x, 0, netlist[u].c, netlist[u].d, hPrint)) return ERRO;
    }
    else if (tipo == 'H') {
      nv = nv + 2;
      if (testarNos(hPrint)) return ERRO;
      strcpy(lista[nv - 1], "jx"); strcat(lista[nv - 1], netlist[u].nome);
      netlist[u].x = nv - 1;
      strcpy(lista[nv], "jy"); strcat(lista[nv], netlist[u].nome);
      netlist[u].y = nv;
      if (operacional(netlist[u].a, netlist[u].b, 0, netlist[u].y, hPrint) ||
        operacional(netlist[u].x, 0, netlist[u].c, netlist[u].d, hPrint)) return ERRO;
    }
    else if (tipo == 'O') {
      if (operacional(netlist[u].a, netlist[u].b, netlist[u].c, netlist[u].d, hPrint)) return ERRO;
      neq--;
    }
    else if (tipo == 'K') {
      nv++;
      neq++; // Como queremos calcular a corrente, nao usamos ampops
      if (testarNos(hPrint)) return ERRO;
      strcpy(lista[nv], "j");
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
    }
    else if (tipo == 'L') {
      nv++;
      neq++;
      if (testarNos(hPrint)) return ERRO;
      strcpy(lista[nv], "j");
      strcat(lista[nv], netlist[u].nome);
      netlist[u].x = nv;
    }
  }

  return 0;
}


/* Rotina que lista as variaveis e o netlist interno final. */
void
listarTudo(HWND hPrint) {
  my_printf(hPrint, "\r\nVariáveis internas:\r\n");
  for (u = 0; u <= nv; u++) {
    sprintf(aux, "  - %d -> %s (%d)\r\n", u, lista[u], C[u]);
    my_printf(hPrint, aux);
  }

  my_printf(hPrint, "\r\nNetlist interno final:\r\n");
  for (u = 1; u <= ne; u++) {
    tipo = netlist[u].nome[0];
    if (tipo == 'R' || tipo == 'I' || tipo == 'V' || tipo == 'C' || tipo == 'L') {
      sprintf(aux, "  - %s %d %d %g\r\n", netlist[u].nome, netlist[u].a, netlist[u].b, netlist[u].valor);
      my_printf(hPrint, aux);
    }
    else if (tipo == 'G' || tipo == 'E' || tipo == 'F' || tipo == 'H') {
      sprintf(aux, "  - %s %d %d %d %d %g\r\n", netlist[u].nome, netlist[u].a, netlist[u].b,
        netlist[u].c, netlist[u].d, netlist[u].valor);
      my_printf(hPrint, aux);
    }
    else if (tipo == 'O') {
      sprintf(aux, "  - %s %d %d %d %d\r\n", netlist[u].nome, netlist[u].a, netlist[u].b,
        netlist[u].c, netlist[u].d);
      my_printf(hPrint, aux);
    }

    if (tipo == 'V' || tipo == 'E' || tipo == 'F' || tipo == 'O' || tipo == 'K' || tipo == 'L') {
      sprintf(aux, "  - Corrente jx [não calculada]: %d\r\n", netlist[u].x);
      my_printf(hPrint, aux);
    }
    else if (tipo == 'H') {
      sprintf(aux, "  - Correntes jx e jy [não calculada]: %d, %d\r\n", netlist[u].x, netlist[u].y);
      my_printf(hPrint, aux);
    }
  }
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


/* Funcao que faz a analise. No programa original, sem a interface grafica, essa e sua main(). */
int
main_cmna(HWND hPrint) {
  srand((unsigned int)time(NULL));

  for (u = 0; u <= MAX_NOS; u++) { // Inicializacao de tabelas
    C[u] = u;
    L[u] = u;
    t0[u] = 0.0;
  }

  if (lerNetlist(hPrint)) return ERRO; // Chamada da rotina que le o netlist
  if (elementosModificada(hPrint)) return ERRO; // Processamento de elementos da analise modificada
  listarTudo(hPrint); // Listagem de variaveis e elementos

  passo = passo / passosInt;
  qtdePontos = (int)round(tempo / passo);
  correcaoPulse();

  strcpy(nomeValores, nomeArquivo); // Cria o nome do arquivo de t0, <nomeArquivo>.tab
  char *pExtensao = strrchr(nomeValores, '.');
  strcpy(pExtensao, ".tab");

  // Inicio da analise
  sprintf(aux, "\r\nO circuito tem %d nós, %d variáveis internas, %d equações e %d elementos.\r\n", nn, nv, neq, ne);
  my_printf(hPrint, aux);

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
      zerarSistema();
      montarEstampas();

      /* Esse while diminui a tolerancia caso o sistema seja singular, porque as vezes
      * a singularidade ocorre por conta de erros numericos ou valores ruins no NR. */
      while (resolverSistema()) {
        if (tolg > TOLG_MIN)
          tolg *= 1e-1;

        if (tolg < TOLG_MIN) {
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
      sprintf(aux, "\r\n(!) ERRO: O sistema não convergiu no ponto %d.\r\n", ponto);
      my_printf(hPrint, aux);
      fclose(valores);
      return ERRO;
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
  sprintf(aux, "\r\nPronto. %d pontos calculados internamente; %d foram incluidos na tabela.\r\n",
    ponto - 1, (ponto - 1) / passosInt);
  my_printf(hPrint, aux);
  sprintf(aux, "O programa demorou %.4Lg s para simular o circuito.\r\n",
    (double)(fim - inicio) / CLOCKS_PER_SEC);
  my_printf(hPrint, aux);
  fclose(valores);

  return 0;
}


/* Rotina da janela, que recebe inputs do usuario. */
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
  switch (msg)
  {
  // Caso base de inicalizacao do programa
  case WM_CREATE:
  {
    HFONT hfDefault;
    HWND hEdit;

    /* Criamos a janela de log. Ela sera de edicao, para que o usuario possa copiar as mensagens,
     * mas nao sera editavel. Como ha tratamento para o tamanho da janela, colocamos qualquer coisa. */
    hEdit = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "",
      WS_CHILD | WS_VISIBLE | WS_VSCROLL | ES_MULTILINE | ES_AUTOVSCROLL | ES_READONLY,
      0, 0, 100, 100, hwnd, (HMENU)IDC_MAIN_EDIT, GetModuleHandle(NULL), NULL);
    if (hEdit == NULL)
      MessageBox(hwnd, "Não foi possível criar a caixa de log.", "Erro!", MB_OK | MB_ICONERROR);

    hfDefault = (HFONT)GetStockObject(DEFAULT_GUI_FONT);
    SendMessage(hEdit, WM_SETFONT, (WPARAM)hfDefault, MAKELPARAM(FALSE, 0));

    int idx = GetWindowTextLength(hEdit);
    SendMessage(hEdit, EM_SETSEL, idx, idx);

    str = "cMNA - Análise Nodal Modificada Compacta no Tempo (Versão 1.0 - 2018)\r\n"
      "Por Matheus F. Moreno e Paulo Victor S. Lima\r\n"
      "Código-base por Antônio Carlos M. de Queiroz; ícone cortesia de Instructables\r\n";

    SendMessage(hEdit, EM_REPLACESEL, 0, (LPARAM)str.c_str());
  }
  break;

  // Reconfiguracao da janela caso o usuario mude seu tamanho
  case WM_SIZE:
  {
    HWND hEdit;
    RECT rcClient;

    GetClientRect(hwnd, &rcClient);

    hEdit = GetDlgItem(hwnd, IDC_MAIN_EDIT);
    SetWindowPos(hEdit, NULL, 0, 0, rcClient.right, rcClient.bottom, SWP_NOZORDER);
  }
  break;

  case WM_COMMAND:
  {
    HWND hEdit = GetDlgItem(hwnd, IDC_MAIN_EDIT);
    switch (LOWORD(wParam))
    {
      // Leitura do nome do arquivo. Esse bloco nao abre o arquivo.
      case ID_ARQUIVO_ABRIR:
      {
        OPENFILENAME ofn;

        ZeroMemory(&ofn, sizeof(ofn));

        ofn.lStructSize = sizeof(ofn);
        ofn.hwndOwner = hwnd;
        ofn.lpstrFilter = "Netlists (*.net)\0*.net\0Todos os arquivos (*.*)\0*.*\0";
        ofn.lpstrFile = nomeArquivo;
        ofn.nMaxFile = MAX_PATH;
        ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
        ofn.lpstrDefExt = "net";

        // Tentamos pegar o nome do arquivo e, se conseguirmos, rodamos a main original
        if (GetOpenFileName(&ofn))
          if (main_cmna(hEdit))
            my_printf(hEdit, "Cheque se o circuito está correto e tente novamente.\r\n");
      }
        break;

      // Sair do programa pela opcao Sair no menu
      case ID_ARQUIVO_SAIR:
        PostMessage(hwnd, WM_CLOSE, 0, 0);
        break;

      // Reanalisar o ultimo circuito. O nome padrao e "circuito.net"
      case ID_REANALISAR:
        if (main_cmna(hEdit))
          my_printf(hEdit, "Cheque se o circuito está correto e tente novamente.\r\n");
        break;

      // Escreve na tela o texto de ajuda
      case ID_AJUDA:
      {
        int idx = GetWindowTextLength(hEdit);
        SendMessage(hEdit, EM_SETSEL, idx, idx);

        SendMessage(hEdit, EM_REPLACESEL, 0, (LPARAM)ajuda.c_str());
      }
        break;
    }
  }
    break;

  case WM_CLOSE:
    DestroyWindow(hwnd);
    break;

  case WM_DESTROY:
    PostQuitMessage(0);
    break;

  default:
    return DefWindowProc(hwnd, msg, wParam, lParam);
  }
  return 0;
}


/* Main do programa. Observe que nao e um main(), e sim um WinMain. */
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
  LPSTR lpCmdLine, int nCmdShow)
{
  WNDCLASSEX wc;
  HWND hwnd;
  MSG Msg;

  // Registrando a classe da janela
  wc.cbSize = sizeof(WNDCLASSEX);
  wc.style = 0;
  wc.lpfnWndProc = WndProc;
  wc.cbClsExtra = 0;
  wc.cbWndExtra = 0;
  wc.hInstance = hInstance;
  wc.hCursor = LoadCursor(NULL, IDC_ARROW);
  wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
  wc.lpszMenuName = NULL;
  wc.lpszClassName = g_szClassName;
  wc.lpszMenuName = MAKEINTRESOURCE(IDR_MENU);
  wc.hIcon = LoadIcon(GetModuleHandle(NULL), MAKEINTRESOURCE(IDI_ICONE));
  wc.hIconSm = (HICON)LoadImage(GetModuleHandle(NULL), MAKEINTRESOURCE(IDI_ICONE), IMAGE_ICON, 32, 32, 0);

  if (!RegisterClassEx(&wc))
  {
    MessageBox(NULL, "Falha no Registro da Janela!", "Erro!",
      MB_ICONEXCLAMATION | MB_OK);
    return 0;
  }

  // Criando a janela
  hwnd = CreateWindowEx(
    WS_EX_CLIENTEDGE,
    g_szClassName,
    "cMNA - Análise Nodal Modificada Compacta",
    WS_OVERLAPPEDWINDOW,
    CW_USEDEFAULT, CW_USEDEFAULT, 650, 450,
    NULL, NULL, hInstance, NULL);

  if (hwnd == NULL)
  {
    MessageBox(NULL, "Falha na Criação da Janela!", "Erro!",
      MB_ICONEXCLAMATION | MB_OK);
    return 0;
  }

  ShowWindow(hwnd, nCmdShow);
  UpdateWindow(hwnd);

  // Loop de mensagem
  while (GetMessage(&Msg, NULL, 0, 0) > 0)
  {
    TranslateMessage(&Msg);
    DispatchMessage(&Msg);
  }
  return Msg.wParam;
}
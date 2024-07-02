#pragma once
//Copyright Takuya OOURA, 1996-2001
/*
���\�[�X�F
 * ���s��w�������̑�Y��Ǝ����t���[�\�t�g�Ƃ��Ē񋟂���
 * �u�ėp FFT (���� �t�[���G/�R�T�C��/�T�C�� �ϊ�) �p�b�P�[�W�v
 * (http://www.kurims.kyoto-u.ac.jp/~ooura/fft-j.html)
 * ��fft4g.c��蔲��,�ꕔ���ς��s����
*/


/*
n :�f�[�^�� ��̗ݏ撷�ɂ���
a :�f�[�^�Q

ip:int�^�|�C���^
w :cos/sin�e�[�u��
��ip[0]=0�̂ɏ����������
*/

void cdft(int n, int isgn, double* a, int* ip, double* w);
void makewt(int nw, int* ip, double* w);

/* -------- child routines -------- */

void bitrv2(int n, int* ip, double* a);
void bitrv2conj(int n, int* ip, double* a);
void cftfsub(int n, double* a, double* w);
void cftbsub(int n, double* a, double* w);
void cft1st(int n, double* a, double* w);
void cftmdl(int n, int l, double* a, double* w);




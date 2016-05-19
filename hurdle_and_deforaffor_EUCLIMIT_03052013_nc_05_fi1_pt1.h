// 19.07.2011 Deforestation rate for Romania is returned to the previous value
//(10.11.2010) The same as 24.19 but affor and defor rate coefficients and hurdle are changed because "original" forest area is used (4.8 mil ha)
// The same as 24.17 except Romania. The coefficients for Romania are changed to 
// take into account correct deforestation rate 4800 ha/year (email from Roberto Pilli (JRC) 09 Mar 2010 )
// in the version 24.19 Romania and the Netherlands rates are correcteed
void hurdle_aff_deff(void)
{

for (int i=0;i<=244;i++){
Hurdle_opt[i]=1.5;
}

for (int i=0;i<=244;i++) {
    afforRate_opt[i]=2.*0.66;
    deforRate_opt[i]=3.09*0.83;};   


Hurdle_opt[	0	]=	0	;
Hurdle_opt[	1	]=	22.1579089	;
Hurdle_opt[	2	]=	6.882519531	;
Hurdle_opt[	3	]=	0.001	;
Hurdle_opt[	4	]=	0.001	;
Hurdle_opt[	5	]=	51.89443506	;
Hurdle_opt[	6	]=	2.9694	;
Hurdle_opt[	7	]=	0.001	;
Hurdle_opt[	8	]=	1.8	;
Hurdle_opt[	9	]=	3.440240732	;
Hurdle_opt[	10	]=	1.6101	;
Hurdle_opt[	11	]=	0.001	;
Hurdle_opt[	12	]=	0.001	;
Hurdle_opt[	13	]=	0.001	;
Hurdle_opt[	14	]=	1.086	;
Hurdle_opt[	15	]=	1.114960477	;
Hurdle_opt[	16	]=	0.275	;	//0.2	;	//0.121	;	//0.33783462	;
Hurdle_opt[	17	]=	1.397304927	;
Hurdle_opt[	18	]=	1.3444	;
Hurdle_opt[	19	]=	0.4029	;
Hurdle_opt[	20	]=	0.1	;
Hurdle_opt[	21	]=	2.466767188	;
Hurdle_opt[	22	]=	0.017	;
Hurdle_opt[	23	]=	10	;	//4.8	;	//4.3096 ; //52.38785	;
Hurdle_opt[	24	]=	0.8166	;
Hurdle_opt[	25	]=	50.5537	;
Hurdle_opt[	26	]=	7.454823144	;
Hurdle_opt[	27	]=	51.3107	;
Hurdle_opt[	28	]=	0.2695	;
Hurdle_opt[	29	]=	0.001	;
Hurdle_opt[	30	]=	0.4058	;
Hurdle_opt[	31	]=	2.064971447	;
Hurdle_opt[	32	]=	1.4175	;
Hurdle_opt[	33	]=	0.24275	;
Hurdle_opt[	34	]=	60.92290763	;
Hurdle_opt[	35	]=	0.001	;
Hurdle_opt[	36	]=	0.2	;
Hurdle_opt[	37	]=	11.00634766	;
Hurdle_opt[	38	]=	0.848170756	;
Hurdle_opt[	39	]=	0.001	;
Hurdle_opt[	40	]=	2.4446	;
Hurdle_opt[	41	]=	3.799665933	;
Hurdle_opt[	42	]=	2.930394531	;
Hurdle_opt[	43	]=	10	;
Hurdle_opt[	44	]=	1.9308	;
Hurdle_opt[	45	]=	1.0449	;
Hurdle_opt[	46	]=	3.320829688	;
Hurdle_opt[	47	]=	0.001	;
Hurdle_opt[	48	]=	2.330083105	;
Hurdle_opt[	49	]=	0.9942	;
Hurdle_opt[	50	]=	0.001	;
Hurdle_opt[	51	]=	100.4113511	;
Hurdle_opt[	52	]=	4.875186133	;
Hurdle_opt[	53	]=	0.001	;
Hurdle_opt[	54	]=	0.001	;
Hurdle_opt[	55	]=	50.6906	;
Hurdle_opt[	56	]=	51.33730409	;
Hurdle_opt[	57	]=	0.194737754	;
Hurdle_opt[	58	]=	1.5652	;
Hurdle_opt[	59	]=	0.9593	;
Hurdle_opt[	60	]=	0.6475	;
Hurdle_opt[	61	]=	0.928376731	;
Hurdle_opt[	62	]=	0.6	;
Hurdle_opt[	63	]=	1.48	;
Hurdle_opt[	64	]=	1.8	;
Hurdle_opt[	65	]=	2.1086	;
Hurdle_opt[	66	]=	1.1767	;
Hurdle_opt[	67	]=	52.25	;
Hurdle_opt[	68	]=	10	;	//8	;	//5;//3	;	//2.66	;	//1.3824	;
Hurdle_opt[	69	]=	0.9066	;
Hurdle_opt[	70	]=	0.590626409	;
Hurdle_opt[	71	]=	2.5919	;
Hurdle_opt[	72	]=	1.6582	;
Hurdle_opt[	73	]=	0.8	;	//1	;	//1.5	;	//5	;	//0.388730148	;
Hurdle_opt[	74	]=	1.1448	;
Hurdle_opt[	75	]=	0.001	;
Hurdle_opt[	76	]=	5.742746335	;
Hurdle_opt[	77	]=	0.419068443	;
Hurdle_opt[	78	]=	12.24683312	;
Hurdle_opt[	79	]=	1.702824433	;
Hurdle_opt[	80	]=	0.001	;
Hurdle_opt[	81	]=	5.243797136	;
Hurdle_opt[	82	]=	0.4643	;
Hurdle_opt[	83	]=	3.6	;
Hurdle_opt[	84	]=	0.4742	;
Hurdle_opt[	85	]=	0.3	;
Hurdle_opt[	86	]=	15	;
Hurdle_opt[	87	]=	1.769	;
Hurdle_opt[	88	]=	1.4928	;
Hurdle_opt[	89	]=	2.020729775	;
Hurdle_opt[	90	]=	1.5	;
Hurdle_opt[	91	]=	0.001	;
Hurdle_opt[	92	]=	11.34419004	;
Hurdle_opt[	93	]=	0.001	;
Hurdle_opt[	94	]=	0.001	;
Hurdle_opt[	95	]=	0.0775	;
Hurdle_opt[	96	]=	5.049881101	;
Hurdle_opt[	97	]=	0.3084	;
Hurdle_opt[	98	]=	10	;	//3	;	//51.26248936	;
Hurdle_opt[	99	]=	0.0639	;
Hurdle_opt[	100	]=	55	;
Hurdle_opt[	101	]=	0.001	;
Hurdle_opt[	102	]=	1;	//2	;	//10	;	//50	;	//100	;	//0.547863223	;
Hurdle_opt[	103	]=	4.37680357	;
Hurdle_opt[	104	]=	5.885030361	;
Hurdle_opt[	105	]=	1.8	;
Hurdle_opt[	106	]=	2.223188333	;
Hurdle_opt[	107	]=	3.9654	;
Hurdle_opt[	108	]=	0.7301	;
Hurdle_opt[	109	]=	1.6016	;
Hurdle_opt[	110	]=	0.864073441	;
Hurdle_opt[	111	]=	0.64825	;
Hurdle_opt[	112	]=	3.362659375	;
Hurdle_opt[	113	]=	15	;
Hurdle_opt[	114	]=	0.16	;
Hurdle_opt[	115	]=	0.001	;
Hurdle_opt[	116	]=	1.6259	;
Hurdle_opt[	117	]=	0.4199305	;
Hurdle_opt[	118	]=	1.1017	;
Hurdle_opt[	119	]=	0.7507	;
Hurdle_opt[	120	]=	1.3858	;
Hurdle_opt[	121	]=	1.9003	;
Hurdle_opt[	122	]=	50.5337	;
Hurdle_opt[	123	]=	1.1771	;
Hurdle_opt[	124	]=	2.4695	;
Hurdle_opt[	125	]=	0.011	;
Hurdle_opt[	126	]=	10.29551629	;
Hurdle_opt[	127	]=	4.579038636	;
Hurdle_opt[	128	]=	0.621330752	;
Hurdle_opt[	129	]=	11.82316357	;
Hurdle_opt[	130	]=	0.001	;
Hurdle_opt[	131	]=	4.807571078	;
Hurdle_opt[	132	]=	0.8666	;
Hurdle_opt[	133	]=	50.96025	;
Hurdle_opt[	134	]=	4.17625	;
Hurdle_opt[	135	]=	0.001	;
Hurdle_opt[	136	]=	3.902425551	;
Hurdle_opt[	137	]=	0.001	;
Hurdle_opt[	138	]=	32.27099219	;
Hurdle_opt[	139	]=	1.9	;
Hurdle_opt[	140	]=	1.6184	;
Hurdle_opt[	141	]=	1	;
Hurdle_opt[	142	]=	0.001	;
Hurdle_opt[	143	]=	3.505	;
Hurdle_opt[	144	]=	1.7492	;
Hurdle_opt[	145	]=	6.543631167	;
Hurdle_opt[	146	]=	0.9491	;
Hurdle_opt[	147	]=	0.001	;
Hurdle_opt[	148	]=	1.1162	;
Hurdle_opt[	149	]=	1.1549	;
Hurdle_opt[	150	]=	8.370605469	;
Hurdle_opt[	151	]=	0.235345	;
Hurdle_opt[	152	]=	0.001	;
Hurdle_opt[	153	]=	0.1	;
Hurdle_opt[	154	]=	0.7237	;
Hurdle_opt[	155	]=	0.001	;
Hurdle_opt[	156	]=	0.001	;
Hurdle_opt[	157	]=	0.05	;
Hurdle_opt[	158	]=	0.32	;
Hurdle_opt[	159	]=	0.001	;
Hurdle_opt[	160	]=	0.369	;
Hurdle_opt[	161	]=	2.5	;
Hurdle_opt[	162	]=	0.03	;
Hurdle_opt[	163	]=	0.001	;
Hurdle_opt[	164	]=	0.472711541	;
Hurdle_opt[	165	]=	3.3679	;
Hurdle_opt[	166	]=	0.99918775	;
Hurdle_opt[	167	]=	1.923367741	;
Hurdle_opt[	168	]=	0.001	;
Hurdle_opt[	169	]=	3.129882813	;
Hurdle_opt[	170	]=	3.0161	;
Hurdle_opt[	171	]=	0.001	;
Hurdle_opt[	172	]=	3.258247803	;
Hurdle_opt[	173	]=	2.544741006	;
Hurdle_opt[	174	]=	0.4922	;
Hurdle_opt[	175	]=	0.3432	;
Hurdle_opt[	176	]=	51	;
Hurdle_opt[	177	]=	1.8	;
Hurdle_opt[	178	]=	1.4778	;
Hurdle_opt[	179	]=	0.001	;
Hurdle_opt[	180	]=	1.5635	;
Hurdle_opt[	181	]=	2.9979	;
Hurdle_opt[	182	]=	4.983686523	;
Hurdle_opt[	183	]=	1.902494315	;
Hurdle_opt[	184	]=	2.5	;
Hurdle_opt[	185	]=	4.9443	;
Hurdle_opt[	186	]=	3	;
Hurdle_opt[	187	]=	4.831695505	;
Hurdle_opt[	188	]=	0.418654832	;
Hurdle_opt[	189	]=	0.001	;
Hurdle_opt[	190	]=	0.001	;
Hurdle_opt[	191	]=	1.2501	;
Hurdle_opt[	192	]=	0.90875	;
Hurdle_opt[	193	]=	0.984	;
Hurdle_opt[	194	]=	0.005	;
Hurdle_opt[	195	]=	2.2965	;
Hurdle_opt[	196	]=	0.03	;
Hurdle_opt[	197	]=	0.001	;
Hurdle_opt[	198	]=	52.08009473	;
Hurdle_opt[	199	]=	0.3324	;
Hurdle_opt[	200	]=	0.668201355	;
Hurdle_opt[	201	]=	4.778976662	;
Hurdle_opt[	202	]=	1.269243175	;
Hurdle_opt[	203	]=	1	;	//1.5	;	//2.5	;
Hurdle_opt[	204	]=	50.75	;
Hurdle_opt[	205	]=	0.001	;
Hurdle_opt[	206	]=	4.66748818	;
Hurdle_opt[	207	]=	0.8661	;
Hurdle_opt[	208	]=	0.2242	;
Hurdle_opt[	209	]=	1.1276	;
Hurdle_opt[	210	]=	0.47	;
Hurdle_opt[	211	]=	0.5825	;
Hurdle_opt[	212	]=	0.001	;
Hurdle_opt[	213	]=	1.6436	;
Hurdle_opt[	214	]=	0.001	;
Hurdle_opt[	215	]=	0.001	;
Hurdle_opt[	216	]=	0.3043	;
Hurdle_opt[	217	]=	4.4706	;
Hurdle_opt[	218	]=	8.044921875	;
Hurdle_opt[	219	]=	0.001	;
Hurdle_opt[	220	]=	0.001	;
Hurdle_opt[	221	]=	5.501582686	;
Hurdle_opt[	222	]=	0.6	;
Hurdle_opt[	223	]=	12.04296875	;
Hurdle_opt[	224	]=	0.001	;
Hurdle_opt[	225	]=	3.327457173	;
Hurdle_opt[	226	]=	2.025208984	;
Hurdle_opt[	227	]=	1.4734	;
Hurdle_opt[	228	]=	0.001	;
Hurdle_opt[	229	]=	2.1204	;
Hurdle_opt[	230	]=	0.2018	;
Hurdle_opt[	231	]=	0.001	;
Hurdle_opt[	232	]=	0.001	;
Hurdle_opt[	233	]=	7.910984375	;
Hurdle_opt[	234	]=	1.7931	;
Hurdle_opt[	235	]=	0.001	;
Hurdle_opt[	236	]=	0.001	;
Hurdle_opt[	237	]=	1.3946	;
Hurdle_opt[	238	]=	5.5722548	;
Hurdle_opt[	239	]=	1	;
Hurdle_opt[	240	]=	0.1788	;
Hurdle_opt[	241	]=	0.001	;
Hurdle_opt[	242	]=	0.001	;
Hurdle_opt[	243	]=	0.001	;
				
afforRate_opt[	0	]=	1	;
afforRate_opt[	1	]=	1	;
afforRate_opt[	2	]=	0.001	;
afforRate_opt[	3	]=	1	;
afforRate_opt[	4	]=	1	;
afforRate_opt[	5	]=	0.753127458	;
afforRate_opt[	6	]=	1	;
afforRate_opt[	7	]=	1	;
afforRate_opt[	8	]=	1	;
afforRate_opt[	9	]=	7.539035596	;
afforRate_opt[	10	]=	4.469158463	;
afforRate_opt[	11	]=	1	;
afforRate_opt[	12	]=	1	;
afforRate_opt[	13	]=	1	;
afforRate_opt[	14	]=	1	;
afforRate_opt[	15	]=	300	;
afforRate_opt[	16	]=	12.27869184;	//12.91448151	;	//11.01	;	//3.647688673	;	

afforRate_opt[	17	]=	1	;
afforRate_opt[	18	]=	1	;
afforRate_opt[	19	]=	1.166978406	;
afforRate_opt[	20	]=	0.001	;
afforRate_opt[	21	]=	1.813904492	;
afforRate_opt[	22	]=	2.137245655	;
afforRate_opt[	23	]=	1.012648066	;	//1.081	;	//1.012642976	;
afforRate_opt[	24	]=	1	;
afforRate_opt[	25	]=	1	;
afforRate_opt[	26	]=	1	;
afforRate_opt[	27	]=	1.251956059	;
afforRate_opt[	28	]=	3.871254191	;
afforRate_opt[	29	]=	1	;
afforRate_opt[	30	]=	3.720367538	;
afforRate_opt[	31	]=	1.998264376	;
afforRate_opt[	32	]=	1	;
afforRate_opt[	33	]=	0.1	;
afforRate_opt[	34	]=	27.47724332	;
afforRate_opt[	35	]=	1	;
afforRate_opt[	36	]=	1	;
afforRate_opt[	37	]=	0.001	;
afforRate_opt[	38	]=	2.613086226	;
afforRate_opt[	39	]=	1	;
afforRate_opt[	40	]=	1.382093746	;
afforRate_opt[	41	]=	10.27603108	;
afforRate_opt[	42	]=	13.15150966	;
afforRate_opt[	43	]=	3.869496656	;
afforRate_opt[	44	]=	0.001	;
afforRate_opt[	45	]=	0.001	;
afforRate_opt[	46	]=	0.224556879	;
afforRate_opt[	47	]=	1	;
afforRate_opt[	48	]=	3.230819029	;
afforRate_opt[	49	]=	1	;
afforRate_opt[	50	]=	1	;
afforRate_opt[	51	]=	8.713892253	;
afforRate_opt[	52	]=	5.308560517	;
afforRate_opt[	53	]=	1	;
afforRate_opt[	54	]=	1	;
afforRate_opt[	55	]=	1	;
afforRate_opt[	56	]=	0.10648416	;
afforRate_opt[	57	]=	5.897963179	;	//11.52755271	;	

afforRate_opt[	58	]=	1	;
afforRate_opt[	59	]=	1	;
afforRate_opt[	60	]=	0.425056782	;
afforRate_opt[	61	]=	1	;
afforRate_opt[	62	]=	300	;
afforRate_opt[	63	]=	29.63962426	;
afforRate_opt[	64	]=	1	;
afforRate_opt[	65	]=	55.58990011	;
afforRate_opt[	66	]=	1	;
afforRate_opt[	67	]=	2.057341017	;
afforRate_opt[	68	]=	1.649714059	;	//1.570843673	;	//0.222698958	; 		

afforRate_opt[	69	]=	0.001	;
afforRate_opt[	70	]=	14.63443155	;
afforRate_opt[	71	]=	1	;
afforRate_opt[	72	]=	1	;
afforRate_opt[	73	]=	13.60859285	;	//13.93927182	;	//15.13984706	;	//20.56180148	;
afforRate_opt[	74	]=	1	;
afforRate_opt[	75	]=	1	;
afforRate_opt[	76	]=	0.01351069	;
afforRate_opt[	77	]=	6.084568249	;	//9.796581979	;	
afforRate_opt[	78	]=	4.031893114	; 

afforRate_opt[	79	]=	0.001	;
afforRate_opt[	80	]=	1	;
afforRate_opt[	81	]=	0.001	;
afforRate_opt[	82	]=	1	;
afforRate_opt[	83	]=	1	;
afforRate_opt[	84	]=	1	;
afforRate_opt[	85	]=	1	;
afforRate_opt[	86	]=	0.356934904	;
afforRate_opt[	87	]=	1	;
afforRate_opt[	88	]=	1	;
afforRate_opt[	89	]=	3.962346181	;
afforRate_opt[	90	]=	1	;
afforRate_opt[	91	]=	1	;
afforRate_opt[	92	]=	1	;
afforRate_opt[	93	]=	1	;
afforRate_opt[	94	]=	1	;
afforRate_opt[	95	]=	18.18553937	;
afforRate_opt[	96	]=	0.652195692	;
afforRate_opt[	97	]=	20.00133462	;
afforRate_opt[	98	]=	0.579131166	;
afforRate_opt[	99	]=	2.057965285	;
afforRate_opt[	100	]=	2.527120862	;
afforRate_opt[	101	]=	1	;
afforRate_opt[	102	]=	4.876487581	;	//3.628131667	;
afforRate_opt[	103	]=	1	;
afforRate_opt[	104	]=	3.095888966	;
afforRate_opt[	105	]=	1	;
afforRate_opt[	106	]=	1	;
afforRate_opt[	107	]=	2.747944348	;
afforRate_opt[	108	]=	1	;
afforRate_opt[	109	]=	1	;
afforRate_opt[	110	]=	0.521015036	;
afforRate_opt[	111	]=	300	;
afforRate_opt[	112	]=	1	;
afforRate_opt[	113	]=	59.10602501	;
afforRate_opt[	114	]=	10.01357202	;
afforRate_opt[	115	]=	1	;
afforRate_opt[	116	]=	1	;
afforRate_opt[	117	]=	37.28668496	;
afforRate_opt[	118	]=	1	;
afforRate_opt[	119	]=	4.625050839	;
afforRate_opt[	120	]=	1	;
afforRate_opt[	121	]=	1	;
afforRate_opt[	122	]=	1	;
afforRate_opt[	123	]=	1	;
afforRate_opt[	124	]=	1	;
afforRate_opt[	125	]=	22.56112411	;
afforRate_opt[	126	]=	1	;
afforRate_opt[	127	]=	0.577421578	;	//0.566526508	;	//4.179213025	;		

afforRate_opt[	128	]=	1.034090799	;
afforRate_opt[	129	]=	0.268389755	;
afforRate_opt[	130	]=	1	;
afforRate_opt[	131	]=	6.230170501	;
afforRate_opt[	132	]=	1	;
afforRate_opt[	133	]=	1.73191201	;
afforRate_opt[	134	]=	1.838484149	;
afforRate_opt[	135	]=	1	;
afforRate_opt[	136	]=	13.81963511	;
afforRate_opt[	137	]=	1	;
afforRate_opt[	138	]=	1.742406351	;
afforRate_opt[	139	]=	21.63561486	;
afforRate_opt[	140	]=	1	;
afforRate_opt[	141	]=	5.364409868	;
afforRate_opt[	142	]=	1	;
afforRate_opt[	143	]=	0.001	;
afforRate_opt[	144	]=	1	;
afforRate_opt[	145	]=	0.001	;
afforRate_opt[	146	]=	1	;
afforRate_opt[	147	]=	1	;
afforRate_opt[	148	]=	1	;
afforRate_opt[	149	]=	1	;
afforRate_opt[	150	]=	10.92845991	;
afforRate_opt[	151	]=	300	;
afforRate_opt[	152	]=	1	;
afforRate_opt[	153	]=	1	;
afforRate_opt[	154	]=	1	;
afforRate_opt[	155	]=	1	;
afforRate_opt[	156	]=	1	;
afforRate_opt[	157	]=	0.001	;
afforRate_opt[	158	]=	2.129484572	;
afforRate_opt[	159	]=	1	;
afforRate_opt[	160	]=	1.587029374	;
afforRate_opt[	161	]=	145.8590242	;
afforRate_opt[	162	]=	205.3179245	;
afforRate_opt[	163	]=	1	;
afforRate_opt[	164	]=	1	;
afforRate_opt[	165	]=	1	;
afforRate_opt[	166	]=	0.001	;
afforRate_opt[	167	]=	6.115301648	;
afforRate_opt[	168	]=	1	;
afforRate_opt[	169	]=	300	;
afforRate_opt[	170	]=	3.179749572	;
afforRate_opt[	171	]=	1	;
afforRate_opt[	172	]=	16.13501303	;
afforRate_opt[	173	]=	0.314551638	;	//	4.626710542	;		

afforRate_opt[	174	]=	1	;
afforRate_opt[	175	]=	29.88159482	;
afforRate_opt[	176	]=	21.62580668	;
afforRate_opt[	177	]=	2.277744831	;
afforRate_opt[	178	]=	1	;
afforRate_opt[	179	]=	1	;
afforRate_opt[	180	]=	1	;
afforRate_opt[	181	]=	1	;
afforRate_opt[	182	]=	0.094122131	;	//1	;
afforRate_opt[	183	]=	0.666294894	;
afforRate_opt[	184	]=	17.50872164	;
afforRate_opt[	185	]=	1	;
afforRate_opt[	186	]=	0.001	;
afforRate_opt[	187	]=	3.286917333	;
afforRate_opt[	188	]=	1	;
afforRate_opt[	189	]=	1	;
afforRate_opt[	190	]=	1	;
afforRate_opt[	191	]=	1	;
afforRate_opt[	192	]=	300	;
afforRate_opt[	193	]=	6.88324323	;
afforRate_opt[	194	]=	1	;
afforRate_opt[	195	]=	1	;
afforRate_opt[	196	]=	1	;
afforRate_opt[	197	]=	1	;
afforRate_opt[	198	]=	1.731800384	;
afforRate_opt[	199	]=	1	;
afforRate_opt[	200	]=	0.447423458	;
afforRate_opt[	201	]=	0.119513808		;	//0.301567173	;
afforRate_opt[	202	]=	0.438856432	;	//1.491435671	;
afforRate_opt[	203	]=	 2.438490541	;	//2.422807293	;	//2.395819184	;	//2.344595746	;	//2.202752158	;	//2.302973705	;	//2.175953235	;	//2.202752158	;	
	
afforRate_opt[	204	]=	73.7635714	;
afforRate_opt[	205	]=	1	;
afforRate_opt[	206	]=	5.48507176	;
afforRate_opt[	207	]=	1	;
afforRate_opt[	208	]=	1	;
afforRate_opt[	209	]=	0.001	;
afforRate_opt[	210	]=	0.240812169	;
afforRate_opt[	211	]=	1	;
afforRate_opt[	212	]=	1	;
afforRate_opt[	213	]=	1	;
afforRate_opt[	214	]=	1	;
afforRate_opt[	215	]=	1	;
afforRate_opt[	216	]=	7.423304122	;
afforRate_opt[	217	]=	300	;
afforRate_opt[	218	]=	3.412445481	;
afforRate_opt[	219	]=	1	;
afforRate_opt[	220	]=	1	;
afforRate_opt[	221	]=	0.001	;
afforRate_opt[	222	]=	0.001	;
afforRate_opt[	223	]=	0.0001;//0.502373914	;
afforRate_opt[	224	]=	1	;
afforRate_opt[	225	]=	300	;
afforRate_opt[	226	]=	0.808398567	;
afforRate_opt[	227	]=	1	;
afforRate_opt[	228	]=	1	;
afforRate_opt[	229	]=	1	;
afforRate_opt[	230	]=	8.641028894	;
afforRate_opt[	231	]=	1	;
afforRate_opt[	232	]=	1	;
afforRate_opt[	233	]=	11.30323551	;
afforRate_opt[	234	]=	1	;
afforRate_opt[	235	]=	1	;
afforRate_opt[	236	]=	1	;
afforRate_opt[	237	]=	1	;
afforRate_opt[	238	]=	0.106503329	;
afforRate_opt[	239	]=	0.001	;
afforRate_opt[	240	]=	0.001	;
afforRate_opt[	241	]=	1	;
afforRate_opt[	242	]=	1	;
afforRate_opt[	243	]=	1	;
				
deforRate_opt[	0	]=	1	;
deforRate_opt[	1	]=	1	;
deforRate_opt[	2	]=	126.9869083	;
deforRate_opt[	3	]=	1	;
deforRate_opt[	4	]=	1	;
deforRate_opt[	5	]=	0.753127458	;
deforRate_opt[	6	]=	1	;
deforRate_opt[	7	]=	1	;
deforRate_opt[	8	]=	1	;
deforRate_opt[	9	]=	7.539035596	;
deforRate_opt[	10	]=	26.53178452	;
deforRate_opt[	11	]=	1	;
deforRate_opt[	12	]=	1	;
deforRate_opt[	13	]=	1	;
deforRate_opt[	14	]=	1	;
deforRate_opt[	15	]=	4.099899363	;
deforRate_opt[	16	]=	16.37332388	;	//14.6	;	//37.97987556	;
deforRate_opt[	17	]=	1	;
deforRate_opt[	18	]=	1	;
deforRate_opt[	19	]=	36.71638557	;
deforRate_opt[	20	]=	27.37443173	;
deforRate_opt[	21	]=	1.813904492	;
deforRate_opt[	22	]=	2.137245655	;
deforRate_opt[	23	]=	0.280700863	;	//1.21	;	//1.012642976	;	

deforRate_opt[	24	]=	1	;
deforRate_opt[	25	]=	1	;
deforRate_opt[	26	]=	1.305176458	;
deforRate_opt[	27	]=	1.251956059	;
deforRate_opt[	28	]=	3.871254191	;
deforRate_opt[	29	]=	1	;
deforRate_opt[	30	]=	3.720367538	;
deforRate_opt[	31	]=	1.999343829	;
deforRate_opt[	32	]=	1	;
deforRate_opt[	33	]=	6.491995607	;
deforRate_opt[	34	]=	27.47724332	;
deforRate_opt[	35	]=	1	;
deforRate_opt[	36	]=	1	;
deforRate_opt[	37	]=	4.933278649	;
deforRate_opt[	38	]=	2.62763573	;
deforRate_opt[	39	]=	1	;
deforRate_opt[	40	]=	1.382093746	;
deforRate_opt[	41	]=	32.12208121	;
deforRate_opt[	42	]=	13.55972834	;
deforRate_opt[	43	]=	0.0971829	;
deforRate_opt[	44	]=	4.822168061	;
deforRate_opt[	45	]=	119.4972058	;
deforRate_opt[	46	]=	6.690916895	;
deforRate_opt[	47	]=	1	;
deforRate_opt[	48	]=	2.298037825	;
deforRate_opt[	49	]=	1	;
deforRate_opt[	50	]=	1	;
deforRate_opt[	51	]=	0.001	;
deforRate_opt[	52	]=	5.308560517	;
deforRate_opt[	53	]=	1	;
deforRate_opt[	54	]=	1	;
deforRate_opt[	55	]=	1	;
deforRate_opt[	56	]=	0.10648416	;
deforRate_opt[	57	]=	3.448058515	;	//15.29206032	;
deforRate_opt[	58	]=	1	;
deforRate_opt[	59	]=	1	;
deforRate_opt[	60	]=	0.887537312	;
deforRate_opt[	61	]=	1	;
deforRate_opt[	62	]=	300	;
deforRate_opt[	63	]=	29.63962426	;
deforRate_opt[	64	]=	1	;
deforRate_opt[	65	]=	55.58990011	;
deforRate_opt[	66	]=	1	;
deforRate_opt[	67	]=	1.744469531	;
deforRate_opt[	68	]=	1.259333789	;	//1.220470471	;	//0.222698958	;	
deforRate_opt[	69	]=	26.62918203	;
deforRate_opt[	70	]=	18.86687445	;
deforRate_opt[	71	]=	1	;
deforRate_opt[	72	]=	1	;
deforRate_opt[	73	]=	20.56180148	;
deforRate_opt[	74	]=	1	;
deforRate_opt[	75	]=	1	;
deforRate_opt[	76	]=	300	;
deforRate_opt[	77	]=	2.367221105	;	//2.906403321	;	//9.796581979	;
deforRate_opt[	78	]=	4.031893114	;
deforRate_opt[	79	]=	8.004930769	;
deforRate_opt[	80	]=	1	;
deforRate_opt[	81	]=	15.13270142	;
deforRate_opt[	82	]=	1	;
deforRate_opt[	83	]=	1	;
deforRate_opt[	84	]=	1	;
deforRate_opt[	85	]=	1	;
deforRate_opt[	86	]=	0.357999835	;
deforRate_opt[	87	]=	1	;
deforRate_opt[	88	]=	1	;
deforRate_opt[	89	]=	4.029762349	;
deforRate_opt[	90	]=	1	;
deforRate_opt[	91	]=	1	;
deforRate_opt[	92	]=	1.014150698	;
deforRate_opt[	93	]=	1	;
deforRate_opt[	94	]=	1	;
deforRate_opt[	95	]=	18.18553937	;
deforRate_opt[	96	]=	1.040520708	;
deforRate_opt[	97	]=	20.00133462	;
deforRate_opt[	98	]=	0.579131166	;
deforRate_opt[	99	]=	2.057965285	;
deforRate_opt[	100	]=	2.343221298	;
deforRate_opt[	101	]=	1	;
deforRate_opt[	102	]=	3.628131667	;
deforRate_opt[	103	]=	1	;
deforRate_opt[	104	]=	3.095888966	;
deforRate_opt[	105	]=	1	;
deforRate_opt[	106	]=	1	;
deforRate_opt[	107	]=	2.747944348	;
deforRate_opt[	108	]=	1	;
deforRate_opt[	109	]=	1	;
deforRate_opt[	110	]=	0.521015036	;
deforRate_opt[	111	]=	300	;
deforRate_opt[	112	]=	1	;
deforRate_opt[	113	]=	59.10602501	;
deforRate_opt[	114	]=	10.01357202	;
deforRate_opt[	115	]=	1	;
deforRate_opt[	116	]=	1	;
deforRate_opt[	117	]=	300	;
deforRate_opt[	118	]=	1	;
deforRate_opt[	119	]=	4.625050839	;
deforRate_opt[	120	]=	1	;
deforRate_opt[	121	]=	1	;
deforRate_opt[	122	]=	1	;
deforRate_opt[	123	]=	1	;
deforRate_opt[	124	]=	1	;
deforRate_opt[	125	]=	22.56112411	;
deforRate_opt[	126	]=	1	;
deforRate_opt[	127	]=	0.080369552	;	//0.091049467	;	//4.179213025	;	

deforRate_opt[	128	]=	42.17958425	;
deforRate_opt[	129	]=	2.46023438	;
deforRate_opt[	130	]=	1	;
deforRate_opt[	131	]=	6.230170501	;
deforRate_opt[	132	]=	1	;
deforRate_opt[	133	]=	1.73191201	;
deforRate_opt[	134	]=	31.61320009	;
deforRate_opt[	135	]=	1	;
deforRate_opt[	136	]=	13.81963511	;
deforRate_opt[	137	]=	1	;
deforRate_opt[	138	]=	1.742406351	;
deforRate_opt[	139	]=	21.63561486	;
deforRate_opt[	140	]=	1	;
deforRate_opt[	141	]=	5.364409868	;
deforRate_opt[	142	]=	1	;
deforRate_opt[	143	]=	300	;
deforRate_opt[	144	]=	1	;
deforRate_opt[	145	]=	41.19705691	;
deforRate_opt[	146	]=	1	;
deforRate_opt[	147	]=	1	;
deforRate_opt[	148	]=	1	;
deforRate_opt[	149	]=	1	;
deforRate_opt[	150	]=	10.92845991	;
deforRate_opt[	151	]=	2.339287029	;
deforRate_opt[	152	]=	1	;
deforRate_opt[	153	]=	1	;
deforRate_opt[	154	]=	1	;
deforRate_opt[	155	]=	1	;
deforRate_opt[	156	]=	1	;
deforRate_opt[	157	]=	300	;
deforRate_opt[	158	]=	2.129484572	;
deforRate_opt[	159	]=	1	;
deforRate_opt[	160	]=	144.4759737	;
deforRate_opt[	161	]=	0.001	;
deforRate_opt[	162	]=	205.3179245	;
deforRate_opt[	163	]=	1	;
deforRate_opt[	164	]=	1	;
deforRate_opt[	165	]=	1	;
deforRate_opt[	166	]=	300	;
deforRate_opt[	167	]=	6.115301648	;
deforRate_opt[	168	]=	1	;
deforRate_opt[	169	]=	7.172281923	;
deforRate_opt[	170	]=	3.179749572	;
deforRate_opt[	171	]=	1	;
deforRate_opt[	172	]=	13.81546167	;
deforRate_opt[	173	]=	0.01078547	;	//0.01273488	;	//4.694463793	;
deforRate_opt[	174	]=	1	;
deforRate_opt[	175	]=	29.88159482	;
deforRate_opt[	176	]=	0.92290851	;
deforRate_opt[	177	]=	2.277744831	;
deforRate_opt[	178	]=	1	;
deforRate_opt[	179	]=	1	;
deforRate_opt[	180	]=	1	;
deforRate_opt[	181	]=	1	;
deforRate_opt[	182	]=	0.172059897	;	//0.175335312	;	//1.086830643	;		

deforRate_opt[	183	]=	0.435294313	;
deforRate_opt[	184	]=	17.50872164	;
deforRate_opt[	185	]=	1	;
deforRate_opt[	186	]=	164.7751685	;
deforRate_opt[	187	]=	3.286917333	;
deforRate_opt[	188	]=	4.836344257	;
deforRate_opt[	189	]=	1	;
deforRate_opt[	190	]=	1	;
deforRate_opt[	191	]=	1	;
deforRate_opt[	192	]=	300	;
deforRate_opt[	193	]=	6.88324323	;
deforRate_opt[	194	]=	1	;
deforRate_opt[	195	]=	1	;
deforRate_opt[	196	]=	1	;
deforRate_opt[	197	]=	1	;
deforRate_opt[	198	]=	1.731800384	;
deforRate_opt[	199	]=	1	;
deforRate_opt[	200	]=	0.836005296	;
deforRate_opt[	201	]=	0.079661514	;	//0.312727435	;	

deforRate_opt[	202	]=	0.237108426	;	//2.250844078	;		

deforRate_opt[	203	]=		18.07184825	;	//19.11720668	;	//21.77438004	;	//28.30615241	;	
		
deforRate_opt[	204	]=	0.001	;
deforRate_opt[	205	]=	1	;
deforRate_opt[	206	]=	5.48507176	;
deforRate_opt[	207	]=	1	;
deforRate_opt[	208	]=	1	;
deforRate_opt[	209	]=	79.34557247	;
deforRate_opt[	210	]=	0.240812169	;
deforRate_opt[	211	]=	1	;
deforRate_opt[	212	]=	1	;
deforRate_opt[	213	]=	1	;
deforRate_opt[	214	]=	1	;
deforRate_opt[	215	]=	1	;
deforRate_opt[	216	]=	7.423304122	;
deforRate_opt[	217	]=	0.001	;
deforRate_opt[	218	]=	3.412445481	;
deforRate_opt[	219	]=	1	;
deforRate_opt[	220	]=	1	;
deforRate_opt[	221	]=	300	;
deforRate_opt[	222	]=	48.82938914	;
deforRate_opt[	223	]=	0.01;//0.502373914	;
deforRate_opt[	224	]=	1	;
deforRate_opt[	225	]=	0.001	;
deforRate_opt[	226	]=	0.826928603	;
deforRate_opt[	227	]=	1	;
deforRate_opt[	228	]=	1	;
deforRate_opt[	229	]=	1	;
deforRate_opt[	230	]=	8.641028894	;
deforRate_opt[	231	]=	1	;
deforRate_opt[	232	]=	1	;
deforRate_opt[	233	]=	12.47062856	;
deforRate_opt[	234	]=	1	;
deforRate_opt[	235	]=	1	;
deforRate_opt[	236	]=	1	;
deforRate_opt[	237	]=	1	;
deforRate_opt[	238	]=	0.109109024	;
deforRate_opt[	239	]=	2.660304637	;
deforRate_opt[	240	]=	42.04896785	;
deforRate_opt[	241	]=	1	;
deforRate_opt[	242	]=	1	;
deforRate_opt[	243	]=	1	;

//--------------------------------------------------------

//---- EU + Cr ---------------------------------------

Hurdle_opt[	16	]=  100;	            // 0.432278001	;	//	0.342600026	;
Hurdle_opt[	19	]=	0.35	;	//	0.27	;
Hurdle_opt[	23	]=	10.30357026	;	//	10	;
Hurdle_opt[	55	]=	50.5	;	//	50.6906	;
Hurdle_opt[	56	]=	3.362964696	;	//	1.8	;
Hurdle_opt[	57	]=	0.461730154	;	//	0.410024607	;
Hurdle_opt[	60	]=	0.719518529	;	//	0.410024607	;
Hurdle_opt[	67	]=	8.585797399	;	//	52.25	;
Hurdle_opt[	68	]=	10.05931638	;	//	4	;
Hurdle_opt[	70	]=	0.63	;	//	0.65	;	//	0.65	;
Hurdle_opt[	73	]=	0.980180819	;	//	1.9	;//2
Hurdle_opt[	77	]=	0.883310904	;	//	0.410024607	;
Hurdle_opt[	86	]=	27.33200874	;	//	25;	//15
Hurdle_opt[	96	]=	5.009774251	;	//	5.055844058	;
Hurdle_opt[	98	]=	5.19398836	;	//	5.19398836	;
Hurdle_opt[	102	]=	1.210359412	;	//	2	;
Hurdle_opt[	107	]=	3.71891531	;	//	3.71891531	;
Hurdle_opt[	127	]=	8.148961633	;	//	5.194146729	;
Hurdle_opt[	128	]=	0.45	;	//	0.45	;
Hurdle_opt[	129	]=	10.80585622	;	//	11.66153243	;
Hurdle_opt[	140	]=	1	;	//	1.6184	;
Hurdle_opt[	160	]=	0.301366051	;	//	0.37	;
Hurdle_opt[	173	]=	3.147248546	;	//	3.5	;
Hurdle_opt[	176	]=	3.418530909	;	//	4	;	//	5.21475625	;
Hurdle_opt[	182	]=	4.916569738	;	//	5.21475625	;
Hurdle_opt[	201	]=	5.39333952	;	//	5.21475625	;
Hurdle_opt[	202	]=	2.017236084	;	//	1.011657889	;
Hurdle_opt[	203	]=	0.915683297	;	//	0.976585938	;
							
afforRate_opt[	16	]=	3.427877073	;
afforRate_opt[	19	]=	1.195875777	;
afforRate_opt[	23	]=	0.914573201	;
afforRate_opt[	55	]=	1	;
afforRate_opt[	56	]=	0.147072301	;
afforRate_opt[	57	]=	0.378589099	;
afforRate_opt[	60	]=	0.390499461	;
afforRate_opt[	67	]=	1.926278747	;
afforRate_opt[	68	]=	1.564767662	;
afforRate_opt[	70	]=	35.8	;
afforRate_opt[	73	]=	4.503878217	;
afforRate_opt[	77	]=	0.443596081	;
afforRate_opt[	86	]=	6.374805842	;
afforRate_opt[	96	]=	1.72862831	;
afforRate_opt[	98	]=	0.560729351	;
afforRate_opt[	102	]=	5.705940581	;	//29.93380852	;
afforRate_opt[	107	]=	3.292867583	;
afforRate_opt[	127	]=	0.55057707	;
afforRate_opt[	128	]=	0.93017563	;
afforRate_opt[	129	]=	0.245303668	;
afforRate_opt[	140	]=	1	;
afforRate_opt[	160	]=	16.3901805	;
afforRate_opt[	173	]=	0.594046181	;
afforRate_opt[	176	]=	18.52074473	;
afforRate_opt[	182	]=	0.102887582	;
afforRate_opt[	201	]=	0.107875719	;
afforRate_opt[	202	]=	1.51649086	;
afforRate_opt[	203	]=	4	;	//8.26	;	//11.80782909	;
				
				
deforRate_opt[	16	]=	16.09329923	;
deforRate_opt[	19	]=	31.76851033	;
deforRate_opt[	23	]=	0.251823988	;
deforRate_opt[	55	]=	1	;
deforRate_opt[	56	]=	0.13814677	;
deforRate_opt[	57	]=	3.396128699	;
deforRate_opt[	60	]=	1.678945004	;
deforRate_opt[	67	]=	1.22230075	;
deforRate_opt[	68	]=	0.426653776	;
deforRate_opt[	70	]=	130	;
deforRate_opt[	73	]=	25.68412361	;
deforRate_opt[	77	]=	5.225045957	;
deforRate_opt[	86	]=	0.238055879	;
deforRate_opt[	96	]=	1.491011367	;
deforRate_opt[	98	]=	0.154598918	;
deforRate_opt[	102	]=	4.85876288	;
deforRate_opt[	107	]=	1.114363016	;
deforRate_opt[	127	]=	0.146961636	;
deforRate_opt[	128	]=	4	;
deforRate_opt[	129	]=	0.194531233	;
deforRate_opt[	140	]=	1	;
deforRate_opt[	160	]=	106.8344711	;
deforRate_opt[	173	]=	0.105912985	;
deforRate_opt[	176	]=	35.987427	;
deforRate_opt[	182	]=	0.167504309	;
deforRate_opt[	201	]=	0.122647001	;
deforRate_opt[	202	]=	0.607737368	;
deforRate_opt[	203	]=	38.86583379	;


//------------------------------------------------------
//------------------------------------------------------


}
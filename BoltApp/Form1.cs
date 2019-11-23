using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace BoltApp
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();


            int nh = 3;//行数
            textBox1.Text = nh.ToString();

            int nv = 3;//列数
            textBox2.Text = nv.ToString();

            double sh = 100;//列距
            textBox4.Text = sh.ToString();

            double sv = 100;//行距
            textBox3.Text = sv.ToString();

            double lx = 200;//锚板x方向长度
            textBox5.Text = lx.ToString();

            double ly = 200;//锚板y方向长度
            textBox6.Text = ly.ToString();

            double tb = 20;//锚板厚度
            textBox7.Text = tb.ToString();

            double Dm = 235;//锚板强度等级
            textBox11.Text = Dm.ToString();

            double Dc = 30;//混凝土基材强度等级
            textBox8.Text = Dc.ToString();

            double d = 20;//锚筋直径
            textBox9.Text = d.ToString();

            double Dl = 4.6;//螺栓等级
            textBox10.Text = Dl.ToString();

            double Nz = 0;//z轴拉力
            textBox12.Text = Nz.ToString();

            double Nzc = -100;//z轴压力,有此项时正压力填为0
            textBox13.Text = Nzc.ToString();

            double Vx = 200;//x向剪力
            textBox14.Text = Vx.ToString();

            double Vy = 0;//y向剪力
            textBox15.Text = Vy.ToString();

            double Mx = 0;//绕x轴正向弯矩
            textBox16.Text = Mx.ToString();

            double My = 200;//绕y轴负向弯矩
            textBox17.Text = My.ToString();

            double Tz = 0;//绕z轴扭矩
            textBox18.Text = Tz.ToString();

        }

        private void printS(double value, bool isMove2NewLine = true)
        {
            if(isMove2NewLine)
                ouput.AppendText(string.Format("{0:#0.0}", value) + nLine);
            else
                ouput.AppendText(string.Format("{0:#0.0}", value) + "\t");
        }

        private void printD(double value1, double value2)
        {
            ouput.AppendText(string.Format("{0:#0.0} \t {0:#0.0} \r\n", value1, value2));
        }
        private void lineFeed( )
        {
            ouput.AppendText(nLine);
        }

        private bool checkValueInt(TextBox tb)
        {
            try
            {
                int i = int.Parse(tb.Text); 
            }
            catch (Exception )
            {
                return false;
            }

            return true;
        }

        private int ParseInt(TextBox tb)
        {
            return int.Parse(tb.Text);
        }

        private bool checkValueDouble(TextBox tb)
        {
            try
            {
                double i = double.Parse(tb.Text);
            }
            catch (Exception)
            {
                return false;
            }

            return true;
        }

        private double ParseDouble(TextBox tb)
        {
            return double.Parse(tb.Text);
        }

        const double PI = 3.1415926;
        string nLine = System.Environment.NewLine;
        private void button1_Click(object sender, EventArgs e)
        {
           BoltAlgorithm bolt = new BoltAlgorithm();

            /*---------------------------------------------------------------------基本参数*/

            int nh = 3;//行数
            if (!checkValueInt(textBox1))
            {
                ouput.AppendText("行数 输入格式错误！");
                return;
            }
            nh = ParseInt(textBox1);

            int nv = 3;//列数
            if (!checkValueInt(textBox2))
            {
                ouput.AppendText("列数 输入格式错误！");
                return;
            }
            nv = ParseInt(textBox2);

            int n = nh * nv;//螺栓个数

            double sh = 100;//列距
            if (!checkValueDouble(textBox4))
            {
                ouput.AppendText("列距 输入格式错误！");
                return;
            }
            sh = ParseDouble(textBox4);

            double sv = 100;//行距
            if (!checkValueDouble(textBox3))
            {
                ouput.AppendText("行距 输入格式错误！");
                return;
            }
            sv = ParseDouble(textBox3);

            /*
                double b1=0;//板边距1
                double b2=0;//板边距2
                double b3=0;//板边距3
                double b4=0;//板边距4
                double c1=0;//基材边距1
                double c2=0;//基材边距2
                double c3=0;//基材边距3
                double c4=0;//基材边距4
            */

            double lx = 200;//锚板x方向长度
            if (!checkValueDouble(textBox5))
            {
                ouput.AppendText("锚板x方向长度 输入格式错误！");
                return;
            }
            lx = ParseDouble(textBox5);

            double ly = 200;//锚板y方向长度
            if (!checkValueDouble(textBox6))
            {
                ouput.AppendText("锚板y方向长度 输入格式错误！");
                return;
            }
            ly = ParseDouble(textBox6);

            double tb = 20;//锚板厚度
            if (!checkValueDouble(textBox7))
            {
                ouput.AppendText("锚板厚度 输入格式错误！");
                return;
            }
            tb = ParseDouble(textBox7);

            double sl = lx * ly;//锚板面积

            double Dm = 235;//锚板强度等级
            if (!checkValueDouble(textBox11))
            {
                ouput.AppendText("锚板强度等级 输入格式错误！");
                return;
            }
            Dm = ParseDouble(textBox11);
            /*
            double tc=300;//混凝土基材厚度
                            */
            double Dc = 30;//混凝土基材强度等级
            if (!checkValueDouble(textBox8))
            {
                ouput.AppendText("混凝土基材强度等级 输入格式错误！");
                return;
            }
            Dc = ParseDouble(textBox8);

            double d = 20;//锚筋直径
            if (!checkValueDouble(textBox9))
            {
                ouput.AppendText("锚筋直径 输入格式错误！");
                return;
            }
            d = ParseDouble(textBox9);

            double Dl = 4.6;//螺栓等级
            if (!checkValueDouble(textBox10))
            {
                ouput.AppendText("螺栓等级 输入格式错误！");
                return;
            }
            Dl = ParseDouble(textBox10);

            double Nz = 0;//z轴拉力
            if (!checkValueDouble(textBox12))
            {
                ouput.AppendText("z轴拉力 输入格式错误！");
                return;
            }
            Nz = ParseDouble(textBox12);

            double Nzc = -100;//z轴压力,有此项时正压力填为0
            if (!checkValueDouble(textBox13))
            {
                ouput.AppendText("z轴压力 输入格式错误！");
                return;
            }
            Nzc = ParseDouble(textBox13);

            double Vx = 200;//x向剪力
            if (!checkValueDouble(textBox14))
            {
                ouput.AppendText("x向剪力 输入格式错误！");
                return;
            }
            Vx = ParseDouble(textBox14);


            double Vy = 0;//y向剪力
            if (!checkValueDouble(textBox15))
            {
                ouput.AppendText("y向剪力 输入格式错误！");
                return;
            }
            Vy = ParseDouble(textBox15);


            double Mx = 0;//绕x轴正向弯矩
            if (!checkValueDouble(textBox16))
            {
                ouput.AppendText("绕x轴正向弯矩 输入格式错误！");
                return;
            }
            Mx = ParseDouble(textBox16);

            double My = 200;//绕y轴负向弯矩
            if (!checkValueDouble(textBox17))
            {
                ouput.AppendText("绕y轴负向弯矩 输入格式错误！");
                return;
            }
            My = ParseDouble(textBox17);

            double Tz = 0;//绕z轴扭矩
            if (!checkValueDouble(textBox18))
            {
                ouput.AppendText("绕z轴扭矩 输入格式错误！");
                return;
            }
            Tz = ParseDouble(textBox18);

            double[] rx = new double[nv];//以螺栓群中心为原点，单个螺栓的x坐标
            double[] ry = new double[nh];//以螺栓群中心为原点，单个螺栓的y坐标
            double[] r = new double[nh * nv];
            double[] rxx = new double[nv];//以螺栓群受压一侧最外排螺栓为原点，单个螺栓的x坐标
            double[] ryy = new double[nh];//以螺栓群受压一侧最外排螺栓为原点，单个螺栓的y坐标
                                         //	printf("依次输入行数，行距，列数，列距，以回车结束\n");
                                         //	scanf("%d,%lf,%d,%lf",&nh,&sv,&nv,&sh);
                                         /*-------------------------------------------------------------轴心拉力基本参数*/
            double Nsd = 0;//螺栓拉力设计值
            double k1 = 1.1;//不均匀系数
                            /*-------------------------------------------------------轴心拉力与弯矩基本参数*/
            double Nsdh = 0;//群锚拉力最大螺栓拉力设计值
            double[] Nsdi = new double[nh * nv];//存放每个螺栓的拉力设计值
            double Nsdg = 0;//群锚受拉区总拉力设计值
                            /*-----------------------------------------------------------剪力与弯矩基本参数*/
            double nx = nh * nv;//参与抵抗x方向剪力螺栓个数
            double ny = nh * nv;//参与抵抗y方向剪力螺栓个数
                                //	double ny=nv;//参与抵抗y方向剪力螺栓个数，边缘，部分抵抗剪力。
            double[] Vvsix = new double[nh * nv];
            double[] Vvsiy = new double[nh * nv];
            double[] Vvs = new double[nh * nv];
            double[] Vtsix = new double[nh * nv];
            double[] Vtsiy = new double[nh * nv];
            double[] Vts = new double[nh * nv];
            double[] Vs = new double[nh * nv];
            double Vsdh = 0.0;
            double ljfh = 0.0;//拉剪复合应力
                        /*-------------------------------------------------------------群锚拉剪许用参数*/
            double nb = 2;//单个螺栓受剪面数
            double Vvb;//单个螺栓受剪承载力设计值
            double Vcb;//单个螺栓承压承载力设计值
            double Ntb;//单个螺栓受拉承载力设计值
                       /*-------------------------------------------------------------群锚拉剪许用限制*/
            Vvb = 0.25 * nb * PI * d * d * bolt.Lfvb(Dl) / 1000;
            Vcb = d * 2 * tb * bolt.Lfcb(Dl, Dm) / 1000;//2*tb，对穿螺栓单个螺栓有2个承压孔壁
            Ntb = 0.25 * PI * bolt.De(d, bolt.P(d)) * bolt.De(d, bolt.P(d)) * bolt.Lftb(Dl) / 10000;
            //	printf("%lf\n",De(d,P(d)));
            /*---------------------------------------计算以螺栓群中心为原点，单个螺栓的坐标*/
            int rxi;
            for (rxi = 0; rxi < nv; rxi++)
            {
                rx[rxi] = bolt.Ljx(nv, sh, rxi + 1);
            }
            /*打印x坐标	
                for(rxi = 0; rxi < nv; rxi++)
                {
                    printf("%lf\n",rx[rxi]);
                }
            */
            int ryi = 0;
            for (ryi = 0; ryi < nh; ryi++)
            {
                ry[ryi] = bolt.Ljy(nh, sv, ryi + 1);
            }
            /*打印y坐标
                for(ryi = 0; ryi < nh; ryi++)
                {
                    printf("%lf\n",ry[ryi]);
                }
            */
            ouput.AppendText("坐标 \r\n");
            for (ryi = 0; ryi < nh; ryi++)
            {
                for (rxi = 0; rxi < nv; rxi++)
                {
                    //printf("(%8.1lf,%8.1lf)\t", rx[rxi], ry[ryi]);
                    printD(rx[rxi], ry[ryi]);
                }
                lineFeed();
                //printf("\n");

            }
            /*-------------------------计算以螺栓群受压一侧最外排螺栓为原点，单个螺栓的坐标*/
            int rxxi;
            for (rxxi = 0; rxxi < nv; rxxi++)
            {
                rxx[rxxi] = bolt.Ljxx(nv, sh, rxxi + 1);
            }
            /*打印xx坐标	
                for(rxxi = 0; rxxi < nv; rxxi++)
                {
                    printf("%lf\n",rxx[rxxi]);
                }
            */
            int ryyi = 0;
            for (ryyi = 0; ryyi < nh; ryyi++)
            {
                ryy[ryyi] = bolt.Ljyy(nh, sv, ryyi + 1);
            }
            /*打印yy坐标
                for(ryyi = 0; ryyi < nh; ryyi++)
                {
                    printf("%lf\n",ryy[ryyi]);
                }
            */
            /*--------------------------计算单个螺栓到螺栓群受压一侧最外排螺栓的x距离平方和*/
            double sumrxxrxx = 0;
            for (rxxi = 0; rxxi < nv; rxxi++)
            {
                sumrxxrxx = sumrxxrxx + rxx[rxxi] * rxx[rxxi];
            }
            /*--------------------------计算单个螺栓到螺栓群受压一侧最外排螺栓的y距离平方和*/
            double sumryyryy = 0;
            for (ryyi = 0; ryyi < nh; ryyi++)
            {
                sumryyryy = sumryyryy + ryy[ryyi] * ryy[ryyi];
            }
            /*-----------------------------------------------计算单个螺栓到螺栓群中心的距离*/
            int ri;
            for (ri = 0; ri < nh * nv; ri++)
            {
                for (ryi = 0; ryi < nh; ryi++)
                {
                    for (rxi = 0; rxi < nv; rxi++)
                    {
                        r[ri] = bolt.Lj(rx[rxi], ry[ryi]);
                        ri++;
                    }
                }
            }
            /*	for (ri  = 0; ri  < nh*nv; ri++)
                {
                    printf("%8.1lf\t",r[ri]);
                    if ((ri+1)%nv==0)
                        printf("%\n");
                }*/
            /*----------------------------------------计算单个螺栓到螺栓群中心的x距离平方和*/
            double sumrxrx = 0;
            for (rxi = 0; rxi < nv; rxi++)
            {
                sumrxrx = sumrxrx + rx[rxi] * rx[rxi];
            }
            /*----------------------------------------计算单个螺栓到螺栓群中心的y距离平方和*/
            double sumryry = 0;
            for (ryi = 0; ryi < nh; ryi++)
            {
                sumryry = sumryry + ry[ryi] * ry[ryi];
            }
            /*-----------------------------------------计算单个螺栓到螺栓群中心的距离平方和*/
            double sumrr = 0;
            for (ri = 0; ri < nh * nv; ri++)
            {
                sumrr = sumrr + r[ri] * r[ri];
            }
            /*---------------------------------------------------群锚在轴心拉力作用下的计算*/
            if ((!(Nz <= 0)) && (Vx == 0) && (Vy == 0) && (Mx == 0) && (My == 0) && (Tz == 0))
            {
                Nsd = k1 * Nz / n;

                ouput.AppendText("群锚受轴心拉力作用" + nLine);
                ouput.AppendText(string.Format("{ 0 :N2}", Nsd)+ nLine);

                //printf("群锚受轴心拉力作用\n");
                //printf("%8.2lf\n", Nsd);
            }
            /*---------------------------------------------群锚在轴心拉力和弯矩作用下的计算*/
            //	if ( ((!(Mx==0))||(!(My==0)))&&(Vx==0)&&(Vy==0)&&(Tz==0))
            //  if ( (!(Mx==0))||(!(My==0)))
            else
            {
                double Mnx;
                double Mny;
                double Mn;
                int Nsdii = 0;
                Mnx = (Nz / n) - ((Mx * ry[0]) / sumryry);
                Mny = (Nz / n) - ((My * rx[nv - 1]) / sumrxrx);
                Mn = (Nz / n) - ((Mx * ry[0]) / sumryry) - ((My * rx[nv - 1]) / sumrxrx);
                if (Mn >= 0)
                {
                    Nsdh = (Nz / n) + ((Mx * ry[0]) / sumryry) + ((My * rx[nv - 1]) / sumrxrx);
                    rxi = 0;
                    ryi = 0;
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        Nsdi[Nsdii] = (Nz / n) + ((Mx * ry[ryi]) / sumryry) + ((My * rx[rxi]) / sumrxrx);
                        //					printf("%8.2lf\n",rxi);
                        //					printf("%8.2lf\n",ryi);
                        rxi++;
                        if ((Nsdii + 1) % nv == 0)
                        {
                            ryi++;
                            rxi = 0;
                        }
                    }
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        Nsdg = Nsdg + Nsdi[Nsdii];
                    }
                    /*			
                                    printf("%8.2lf\n",ry[0]);
                                    printf("%8.2lf\n",rx[nv-1]);
                                    printf("%8.2lf\n",sumryry);
                                    printf("%8.2lf\n",sumrxrx);
                    */
                    //printf("群锚受轴心拉力和弯矩共同作用\n");
                    //printf("拉弯作用下各螺栓拉力设计值\n");
                    ouput.AppendText("拉弯作用下各螺栓拉力设计值" + nLine);
                    //printf("全部螺栓受拉\n");
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        //printf("%8.3lf\t", Nsdi[Nsdii]);
                        printS(Nsdi[Nsdii], false);
                        if ((Nsdii + 1) % nv == 0)
                        {
                            lineFeed();
                            //printf("%\n");
                        }
                    }
                    // printf("Nsdh=%.2lf\n", Nsdh);
                    // printf("Nsdg=%.2lf\n", Nsdg);
                    printS(Nsdh);
                    printS(Nsdg);
                }
                else if ((Mn < 0) && (Mny >= 0) & (Mnx < 0))
                {
                    if (My == 0)
                    {
                        Nsdh = (Nz * bolt.LNy(sv, nh) + Mx) * ryy[0] / sumryyryy;
                        ryi = 0;
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            Nsdi[Nsdii] = (Nz * bolt.LNy(sv, nh) + Mx) * ryy[ryi] / sumryyryy; ;
                            if ((Nsdii + 1) % nv == 0)
                            {
                                ryi++;
                            }
                        }
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            if (Nsdi[Nsdii] >= 0)
                            {
                                Nsdg = Nsdg + Nsdi[Nsdii];
                            }
                        }
                        //printf("群锚受轴心拉力和弯矩共同作用\n");
                        //printf("拉弯作用下各螺栓拉力设计值\n");
                        ouput.AppendText("拉弯作用下各螺栓拉力设计值" + nLine);
                        //printf("绕x轴部分螺栓受拉\n");
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            printS(Nsdi[Nsdii], false);
                            // printf("%8.2lf\t", Nsdi[Nsdii]);
                            if ((Nsdii + 1) % nv == 0)
                                ouput.AppendText(nLine);
                            //printf("%\n");
                        }
                        // printf("Nsdh=%.2lf\n", Nsdh);
                        // printf("Nsdg=%.2lf\n", Nsdg);
                        printS(Nsdh);
                        printS(Nsdg);
                    }
                    else
                    {
                        Nsdh = ((My * rx[nv - 1]) / sumrxrx) + (Nz * bolt.LNy(sv, nh) + Mx) * ryy[0] / sumryyryy;
                        rxi = 0;
                        ryyi = 0;
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            Nsdi[Nsdii] = ((My * rx[rxi]) / sumrxrx) + (Nz * bolt.LNy(sv, nh) + Mx) * ryy[ryyi] / sumryyryy;
                            rxi++;
                            if ((Nsdii + 1) % nv == 0)
                            {
                                rxi = 0;
                                ryyi++;
                            }
                        }
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            if (Nsdi[Nsdii] >= 0)
                            {
                                Nsdg = Nsdg + Nsdi[Nsdii];
                            }
                        }
                        //printf("群锚受轴心拉力和弯矩共同作用\n");
                        // printf("拉弯作用下各螺栓拉力设计值\n");
                        ouput.AppendText("拉弯作用下各螺栓拉力设计值" + nLine);
                        //printf("绕x轴部分螺栓受拉\n");
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            //printf("%8.2lf\t", Nsdi[Nsdii]);
                            printS(Nsdi[Nsdii], false);
                            if ((Nsdii + 1) % nv == 0)
                            {
                                ouput.AppendText(nLine);
                                //printf("%\n");
                            }
                        }
                        //printf("Nsdh=%.2lf\n", Nsdh);
                        //printf("Nsdg=%.2lf\n", Nsdg);
                        printS(Nsdh);
                        printS(Nsdg);
                    }
                }
                else if ((Mn < 0) && (Mnx >= 0) && (Mny < 0))
                {
                    if (Mx == 0)
                    {
                        Nsdh = (Nz * bolt.LNx(sh, nv) + My) * rxx[0] / sumrxxrxx;
                        rxxi = 0;
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            Nsdi[Nsdii] = (Nz * bolt.LNx(sh, nv) + My) * rxx[nv - 1 - rxxi] / sumrxxrxx;
                            rxxi++;
                            if ((Nsdii + 1) % nv == 0)
                            {
                                rxxi = 0;
                            }
                        }
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            if (Nsdi[Nsdii] >= 0)
                            {
                                Nsdg = Nsdg + Nsdi[Nsdii];
                            }
                        }
                        //printf("群锚受轴心拉力和弯矩共同作用\n");
                        //printf("拉弯作用下各螺栓拉力设计值\n");
                        ouput.AppendText("拉弯作用下各螺栓拉力设计值" + nLine);
                        //printf("绕y轴部分螺栓受拉\n");
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            // printf("%8.2lf\t", Nsdi[Nsdii]);
                            printS(Nsdi[Nsdii], false);
                            if ((Nsdii + 1) % nv == 0)
                                ouput.AppendText(nLine);
                            //printf("%\n");
                        }
                        // printf("Nsdh=%.2lf\n", Nsdh);
                        // printf("Nsdg=%.2lf\n", Nsdg);
                        printS(Nsdh);
                        printS(Nsdg);
                    }
                    else
                    {
                        Nsdh = ((Mx * ry[0]) / sumryry) + (Nz * bolt.LNx(sh, nv) + My) * rxx[0] / sumrxxrxx;
                        ryi = 0;
                        rxxi = 0;
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            Nsdi[Nsdii] = ((Mx * ry[ryi]) / sumryry) + (Nz * bolt.LNx(sh, nv) + My) * rxx[nv - 1 - rxxi] / sumrxxrxx;
                            rxxi++;
                            if ((Nsdii + 1) % nv == 0)
                            {
                                rxxi = 0;
                                ryi++;
                            }
                        }
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            if (Nsdi[Nsdii] >= 0)
                            {
                                Nsdg = Nsdg + Nsdi[Nsdii];
                            }
                        }
                        //printf("群锚受轴心拉力和弯矩共同作用\n");
                        //printf("拉弯作用下各螺栓拉力设计值\n");
                        ouput.AppendText("拉弯作用下各螺栓拉力设计值" + nLine);
                        //printf("绕y轴部分螺栓受拉\n");
                        for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                        {
                            // printf("%8.2lf\t", Nsdi[Nsdii]);
                            printS(Nsdi[Nsdii], false);
                            if ((Nsdii + 1) % nv == 0)
                                ouput.AppendText(nLine);
                           // printf("%\n");
                        }
                        // printf("Nsdh=%.2lf\n", Nsdh);
                        // printf("Nsdg=%.2lf\n", Nsdg);
                        printS(Nsdh);
                        printS(Nsdg);
                    }
                }
                else if ((Mn < 0) && (Mnx > 0) && (Mny > 0))
                {
                    Nsdh = (Nz / n) + ((Mx * ry[0]) / sumryry) + ((My * rx[nv - 1]) / sumrxrx);
                    rxi = 0;
                    ryi = 0;
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        Nsdi[Nsdii] = (Nz / n) + ((Mx * ry[ryi]) / sumryry) + ((My * rx[rxi]) / sumrxrx);
                        rxi++;
                        if ((Nsdii + 1) % nv == 0)
                        {
                            ryi++;
                            rxi = 0;
                        }
                    }
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        if (Nsdi[Nsdii] >= 0)
                        {
                            Nsdg = Nsdg + Nsdi[Nsdii];
                        }
                    }
                    //printf("群锚受轴心拉力和弯矩共同作用\n");
                    //printf("拉弯作用下各螺栓拉力设计值\n");
                    ouput.AppendText("拉弯作用下各螺栓拉力设计值" + nLine);
                    //printf("部分螺栓受拉\n");
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        // printf("%8.3lf\t", Nsdi[Nsdii]);  
                        printS(Nsdi[Nsdii], false);
                        if ((Nsdii + 1) % nv == 0)
                            ouput.AppendText(nLine);
                        //printf("%\n");
                    }
                    // printf("Nsdh=%.2lf\n", Nsdh);
                    // printf("Nsdg=%.2lf\n", Nsdg);
                    printS(Nsdh);
                    printS(Nsdg);
                }
                else
                {
                    Nsdh = (Nz * bolt.LNy(sv, nh) + Mx) * ryy[0] / sumryyryy + (Nz * bolt.LNx(sh, nv) + My) * rxx[0] / sumrxxrxx - Nz / n;
                    ryyi = 0;
                    rxxi = 0;
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        Nsdi[Nsdii] = (Nz * bolt.LNy(sv, nh) + Mx) * ryy[ryyi] / sumryyryy + (Nz * bolt.LNx(sh, nv) + My) * rxx[nv - 1 - rxxi] / sumrxxrxx - Nz / n;
                        rxxi++;
                        if ((Nsdii + 1) % nv == 0)
                        {
                            rxxi = 0;
                            ryyi++;
                        }
                    }
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        if (Nsdi[Nsdii] >= 0)
                        {
                            Nsdg = Nsdg + Nsdi[Nsdii];
                        }
                    }
                    //printf("群锚受轴心拉力和弯矩共同作用\n");
                    //printf("拉弯作用下各螺栓拉力设计值\n");
                    ouput.AppendText("拉弯作用下各螺栓拉力设计值" + nLine);
                    //printf("双向部分螺栓受拉\n");
                    for (Nsdii = 0; Nsdii < nh * nv; Nsdii++)
                    {
                        // printf("%8.2lf\t", Nsdi[Nsdii]);
                        printS(Nsdi[Nsdii], false);
                        if ((Nsdii + 1) % nv == 0)
                            ouput.AppendText(nLine);
                        //printf("%\n");
                    }
                    // printf("Nsdh=%.2lf\n", Nsdh);
                    //  printf("Nsdg=%.2lf\n", Nsdg);
                    printS(Nsdh);
                    printS(Nsdg);
                }
            }
            /*-------------------------------------------------群锚在剪力和扭矩作用下的计算*/
            if ((!(Vx == 0)) || (!(Vy == 0)) || (!(Tz == 0)))
            {
                int Vvsi;
                if (nx == ny)
                {
                    for (Vvsi = 0; Vvsi < nh * nv; Vvsi++)
                    {
                        Vvsix[Vvsi] = Vx / nx;
                        Vvsiy[Vvsi] = Vy / ny;
                    }
                }
                else
                {
                    for (Vvsi = 0; Vvsi < nh * nv; Vvsi++)
                    {
                        Vvsix[Vvsi] = Vx / nx;
                        if (Vvsi < ((nh - 1) * nv))
                        {
                            Vvsiy[Vvsi] = 0;
                        }
                        else
                        {
                            Vvsiy[Vvsi] = Vy / ny;
                        }
                    }
                }
                for (Vvsi = 0; Vvsi < nh * nv; Vvsi++)
                {
                    Vvs[Vvsi] = Math.Sqrt((Vvsix[Vvsi]) * (Vvsix[Vvsi]) + (Vvsiy[Vvsi]) * (Vvsiy[Vvsi]));
                }
                int Vtsi = 0;
                /*int*/ rxi = 0;
                /*int*/ ryi = 0;
                for (Vtsi = 0; Vtsi < nh * nv; Vtsi++)
                {
                    Vtsix[Vtsi] = Tz * ry[ryi] / sumrr;
                    Vtsiy[Vtsi] = Tz * rx[rxi] / sumrr;
                    rxi++;
                    if ((Vtsi + 1) % nv == 0)
                    {
                        rxi = 0;
                        ryi++;
                    }
                }
                for (Vtsi = 0; Vtsi < nh * nv; Vtsi++)
                {
                    Vts[Vtsi] = Math.Sqrt((Vtsix[Vtsi]) * (Vtsix[Vtsi]) + (Vtsiy[Vtsi]) * (Vtsiy[Vtsi]));
                }
                int Vsi;
                for (Vsi = 0; Vsi < nh * nv; Vsi++)
                {
                    Vs[Vsi] = Math.Sqrt((Vvsix[Vsi] + Vtsix[Vsi]) * (Vvsix[Vsi] + Vtsix[Vsi]) + (Vvsiy[Vsi] + Vtsiy[Vsi]) * (Vvsiy[Vsi] + Vtsiy[Vsi]));
                }
                for (Vsi = 0; Vsi < nh * nv; Vsi++)
                {
                    if (Vsi == 0)
                    {
                        Vsdh = Vs[Vsi];
                    }
                    else
                        if (Vsdh - Vs[Vsi] < 0)
                    {
                        Vsdh = Vs[Vsi];
                    }
                }
                //printf("群锚受剪力和扭矩共同作用\n");
                /*				printf("剪力作用下各螺栓剪力设计值\n");
                                for (Vvsi = 0; Vvsi < nh*nv; Vvsi++)
                                {
                                    printf("%8.2lf\t",Vvs[Vvsi]);
                                    if ((Vvsi+1)%nv==0)
                                        printf("%\n");
                                }
                                printf("扭矩作用下各螺栓剪力设计值\n");
                                for (Vtsi = 0; Vtsi < nh*nv; Vtsi++)
                                {
                                    printf("%8.2lf\t",Vts[Vtsi]);
                                    if ((Vtsi+1)%nv==0)
                                        printf("%\n");
                                }
                */
                ouput.AppendText("剪扭作用下各螺栓剪力设计值" + nLine);
                //printf("剪扭作用下各螺栓剪力设计值\n");
                for (Vsi = 0; Vsi < nh * nv; Vsi++)
                {
                    // printf("%8.2lf\t", Vs[Vsi]);
                    printS(Vs[Vsi], false);
                    if ((Vsi + 1) % nv == 0)
                        ouput.AppendText(nLine);
                    //printf("%\n");
                }
                // printf("Vsdh=%.2lf\n", Vsdh);
                printS(Vsdh);
                lineFeed();
                /*-----------------------------------------------------------------拉剪复合应力*/
                ljfh = Math.Sqrt((Nsdh / Ntb) * (Nsdh / Ntb) + (Vsdh / Vvb) * (Vsdh / Vvb));
                if (ljfh <= 1 && (Vsdh <= Vcb))
                {
                    // printf("拉剪复合应力=%.2lf<=1.0\n单个螺栓剪力最大值Vsdh=%.2lf<=螺栓承压强度设计值Vcb=%.2lf\n螺栓强度满足要求\n", ljfh, Vsdh, Vcb);
                    ouput.AppendText(string.Format("拉剪复合应力={0:#0.0}<=1.0\r\n单个螺栓剪力最大值Vsdh={0:#0.0} <= 螺栓承压强度设计值Vcb={0:#0.0}\r\n螺栓强度满足要求\r\n", ljfh, Vsdh, Vcb));
                }
                else if (ljfh > 1)
                {
                    // printf("拉剪复合应力=%.2lf>1.0,螺栓强度不满足要求\n", ljfh);
                    ouput.AppendText(string.Format("拉剪复合应力={0:#0.0}>1.0, 螺栓强度不满足要求\r\n", ljfh));
                }
                else if (Vsdh > Vcb)
                {
                    // printf("单个螺栓剪力最大值Vsdh=%.2lf>螺栓承压强度设计值Vcb=%.2lf，螺栓强度不满足要求\n", Vsdh, Vcb);
                    ouput.AppendText(string.Format("单个螺栓剪力最大值Vsdh={0:#0.0} > 螺栓承压强度设计值Vcb={0:#0.0}, 螺栓强度不满足要求\r\n", Vsdh, Vcb));
                }
            }
            /*-----------------------------------------------------------------拉剪复合应力*/
            if (Nzc < 0)
            {
                double fp;
                fp = Math.Abs(Nzc) * 1000 / sl;
                if (fp <= (bolt.Fc(Dc)))
                {
                    //  printf("混凝土在轴压力作用下的强度值fp=%.2lf<=fc=%.2lf,锚板面积满足要求\n", fp, Fc(Dc));
                    ouput.AppendText(string.Format("混凝土在轴压力作用下的强度值 fp={0:#0.0} <= fc={0:#0.0}, 锚板面积满足要求\r\n", fp, bolt.Fc(Dc)));
                }
                else
                {
                    
                    //  printf("混凝土在轴压力作用下的强度值fp=%.2lf>fc=%.2lf,锚板面积不满足要求\n", fp, Fc(Dc));
                    ouput.AppendText(string.Format("混凝土在轴压力作用下的强度值fp={0:#0.0} > fc={0:#0.0}, 锚板面积不满足要求\r\n", fp, bolt.Fc(Dc)));
                }
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            ouput.Clear();
        }
    }
}

vedio = VideoWriter('shock_tube.avi'); %��ʼ��һ��avi�ļ�
vedio.FrameRate = 40;
open(vedio);
for i=1:500  %ͼ�����и���
    fname = [num2str(i,'frame_%05d'),'.jpg'];
    frame = imread(fname); frame = imresize(frame, [1080,1960]);
    writeVideo(vedio,frame);
end
close(vedio);
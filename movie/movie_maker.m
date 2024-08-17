vedio = VideoWriter('shock_tube.avi'); %初始化一个avi文件
vedio.FrameRate = 40;
open(vedio);
for i=1:500  %图像序列个数
    fname = [num2str(i,'frame_%05d'),'.jpg'];
    frame = imread(fname); frame = imresize(frame, [1080,1960]);
    writeVideo(vedio,frame);
end
close(vedio);
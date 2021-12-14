ind1 = 0;
ind2 = 25;


for i = ind2 :-1:ind1
     img = imread(['depth_' num2str(i,'%07d') '.png']);
     filename = ['depth_' num2str(i+1,'%07d') '.png'];
     imwrite(img,filename);
     i
      delete(['depth_' num2str(i,'%07d') '.png'])
end
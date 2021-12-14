ind1 = 0;
ind2 = 25;


for i = ind1 :ind2
     img = imread(['image_' num2str(i,'%07d') '.jpg']);
     filename = [num2str(i+1) '.png'];
     imwrite(img,filename);
     i
     delete(['image_' num2str(i,'%07d') '.jpg'])
end
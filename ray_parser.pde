///////////////////////////////////////////////////////////////////////
//
// Command Line Interface (CLI) Parser  
//
///////////////////////////////////////////////////////////////////////
String gCurrentFile = new String("rect_test.cli"); // A global variable for holding current active file name.

///////////////////////////////////////////////////////////////////////
//
// Press key 1 to 9 and 0 to run different test cases.
//
///////////////////////////////////////////////////////////////////////
void keyPressed() {
  switch(key) {
    case '1':  gCurrentFile = new String("t0.cli"); interpreter(); break;
    case '2':  gCurrentFile = new String("t1.cli"); interpreter(); break;
    case '3':  gCurrentFile = new String("t2.cli"); interpreter(); break;
    case '4':  gCurrentFile = new String("t3.cli"); interpreter(); break;
    case '5':  gCurrentFile = new String("c0.cli"); interpreter(); break;
    case '6':  gCurrentFile = new String("c1.cli"); interpreter(); break;
    case '7':  gCurrentFile = new String("c2.cli"); interpreter(); break;
    case '8':  gCurrentFile = new String("c3.cli"); interpreter(); break;
    case '9':  gCurrentFile = new String("c4.cli"); interpreter(); break;
    case '0':  gCurrentFile = new String("c5.cli"); interpreter(); break;
  }
}

///////////////////////////////////////////////////////////////////////
//
//  Parser core. It parses the CLI file and processes it based on each 
//  token. Only "color", "rect", and "write" tokens are implemented. 
//  You should start from here and add more functionalities for your
//  ray tracer.
//
//  Note: Function "splitToken()" is only available in processing 1.25 
//       or higher.
//
///////////////////////////////////////////////////////////////////////

ArrayList<Renderable> renderables = new ArrayList<Renderable>();
ArrayList<PVector> vertexList = new ArrayList<PVector>();
ArrayList<Light> lights = new ArrayList<Light>();
PVector[][] colorArray = new PVector[300][300];
HitInfo[][] hits = new HitInfo[300][300];

Material currentSurface;
PVector bgColor = new PVector(0,0,0);
int fovDegrees;
float cameraTop;
float cameraBottom;
float cameraLeft;
float cameraRight;
float fovRadians;

int debugCount = 0;

void interpreter() {
  
  String str[] = loadStrings(gCurrentFile);
  if (str == null) println("Error! Failed to read the file.");
  for (int i=0; i<str.length; i++) {
    
    String[] token = splitTokens(str[i], " "); // Get a line and parse tokens.
    if (token.length == 0) continue; // Skip blank line.
    
    if (token[0].equals("fov")) {
      fovDegrees = int(token[1]);
      fovRadians = ((fovDegrees * PI)/180.0);
      cameraTop = (tan(fovRadians/2));
      cameraBottom = -(tan(fovRadians/2));
      cameraLeft = -(tan(fovRadians/2));
      cameraRight = (tan(fovRadians/2));
//      println(cameraTop);
//      println(cameraBottom);
//      println(cameraLeft);
//      println(cameraRight);
    }
    else if (token[0].equals("background")) {
      bgColor = new PVector(float(token[1]), float(token[2]), float(token[3]));
    }
    else if (token[0].equals("light")) {
      PVector position = new PVector(float(token[1]), float(token[2]), float(token[3]));
      PVector intensity = new PVector(float(token[4]), float(token[5]), float(token[6]));
      Light newLight = new Light(position, intensity);
      lights.add(newLight);
    }
    else if (token[0].equals("surface")) {
      PVector diffuse = new PVector(float(token[1]), float(token[2]), float(token[3]));
      PVector ambient = new PVector(float(token[4]), float(token[5]), float(token[6]));
      PVector specular = new PVector(float(token[7]), float(token[8]), float(token[9]));
      int specularExponent = int(token[10]);
      float reflectanceCoefficient = float(token[11]);
      //println(reflectanceCoefficient);
      
      currentSurface = new Material(diffuse, ambient, specular, specularExponent, reflectanceCoefficient);
    }    
    else if (token[0].equals("sphere")) {
      float radius = float(token[1]);
      PVector center = new PVector(float(token[2]), float(token[3]), float(token[4]));
      Material clonedSurface = new Material(currentSurface);
      
      Sphere sphere = new Sphere(radius, center, clonedSurface);
      renderables.add(sphere);
      
    }
    else if (token[0].equals("begin")) {
      // TODO
    }
    else if (token[0].equals("vertex")) {
      PVector vertex = new PVector(float(token[1]), float(token[2]), float(token[3]));
      vertexList.add(vertex);
      if(vertexList.size() == 3)
      {
        Material clonedSurface = new Material(currentSurface);
        Triangle triangle = new Triangle(vertexList.get(0), vertexList.get(1), vertexList.get(2), clonedSurface);
        
        renderables.add((Renderable)triangle);
          //println("We have 3 vertexs");
//          println(vertexList.get(0));
//          println(vertexList.get(1));
//          println(vertexList.get(2));
        vertexList = new ArrayList<PVector>();
      }
    }
    else if (token[0].equals("end")) {
      // TODO
    }
    else if (token[0].equals("color")) {
      float r =float(token[1]);
      float g =float(token[2]);
      float b =float(token[3]);
      fill(r, g, b);
    }
    else if (token[0].equals("rect")) {
      float x0 = float(token[1]);
      float y0 = float(token[2]);
      float x1 = float(token[3]);
      float y1 = float(token[4]);
      rect(x0, height-y1, x1-x0, y1-y0);
    }
    else if (token[0].equals("write")) {
      // you should render the scene here
      
      for(int u = 0; u < height; u++){
        for(int v = 0; v < width; v++){
          
          int actualU = 299 - u;
          
          colorArray[actualU][v] = new PVector(-1, -1, -1);
          
          PVector origin = new PVector(0,0,0);
          PVector target = new PVector(cameraLeft + ((cameraRight - cameraLeft)*(v + .5)/width), cameraBottom + ((cameraTop - cameraBottom)*(u + .5)/height), -1.0);
          PVector direction = new PVector(target.x - origin.x, target.y - origin.y, -1.0);
          PVector normalizedDirection = new PVector(direction.x, direction.y, direction.z);
          normalizedDirection.normalize();          
          
          Ray curRay = new Ray(origin, normalizedDirection);
          
          if(renderables.size() > 0)
          {
            
            int closestRenderable = closestIntersection(renderables, curRay);
            
            if(u == 161 && v == 289){
              println(closestRenderable);
            }
            
            if(closestRenderable == -1 ){
              colorArray[actualU][v] = bgColor;
            }else{
              
              Renderable currentRenderable = renderables.get(closestRenderable);
              HitInfo hitInfo = currentRenderable.intersection(curRay);
              
              hits[actualU][v] = hitInfo;
              
              if(hitInfo.hit)
              {
                
                normalizedDirection.mult(hitInfo.parameterT);
                PVector closestIntersection = normalizedDirection;
                
//                if(hitInfo.triangle)
//                  colorArray[299 - u][v] = computeColor(lights, blockingRenderables, (Triangle)currentRenderable, hitInfo);
//                else
//                  colorArray[299 - u][v] = computeColor(lights, blockingRenderables, (Sphere)currentRenderable, hitInfo, 0);
                colorArray[actualU][v] = computeColor(lights, renderables, currentRenderable, hitInfo, 0);
                
              }else{
                colorArray[actualU][v] = bgColor;
              }
            }
          }
          
        }
      }
      loadPixels();
      if(!gCurrentFile.equals("rect_test.cli")){
        
        for(int u = 0; u < height; u++)
        {
            for (int v = 0; v < width; v++){
              //println(intersectionA1rray[u][v])
              if(colorArray[u][v].x >= 0)
                pixels[u * 300 + v] = color(colorArray[u][v].x, colorArray[u][v].y, colorArray[u][v].z );
              else
                pixels[u * 300 + v] = color(bgColor.x, bgColor.y, bgColor.z );
            }
        }
        
      }
      updatePixels();

      save(token[1]);
      
      renderables.clear();
      lights.clear();
      debugCount = 0;
      bgColor = new PVector(0,0,0);
      
    }
  }
}


int closestIntersection(ArrayList<Renderable> renderables, Ray curRay){
  ArrayList<Integer> intersections = new ArrayList<Integer>();
  for(int k = 0; k < renderables.size(); k++)
  {
    Renderable newRenderable = renderables.get(k);
    HitInfo newHitInfo = newRenderable.intersection(curRay);


    if(newHitInfo.hit)
    {
      intersections.add(k);
    }
  }
  
  debugCount++;
  if(debugCount == (161 * 300) + 289){
    println(intersections.size());
  } 
  
  if(intersections.size() > 0){
    
    Renderable closestRenderable = renderables.get(intersections.get(0));
    HitInfo hitInfo = closestRenderable.intersection(curRay);
    int indexToReturn = intersections.get(0);
    
    for(int k = 1; k < intersections.size(); k++){
      Renderable newRenderable =  renderables.get(intersections.get(k));
      HitInfo newHitInfo = newRenderable.intersection(curRay);
      
      if(newHitInfo.parameterT < hitInfo.parameterT)
      {
        
        closestRenderable = newRenderable;
        hitInfo = newHitInfo;
        indexToReturn = intersections.get(k);
      }
    }
    return indexToReturn;
  }else{
    return -1;
  }
}

PVector computeColor(ArrayList<Light> lights, ArrayList<Renderable> renderables, Renderable curRenderable, HitInfo hitInfo, int recursionLevel){
  //debugCount++;
  if(recursionLevel == 20){
    return new PVector(curRenderable.material.ambient.x, curRenderable.material.ambient.y, curRenderable.material.ambient.z);
  }
  
  //Ambient
  PVector clr = new PVector(curRenderable.material.ambient.x, curRenderable.material.ambient.y, curRenderable.material.ambient.z);

  boolean shadows = false;
  for(int i = 0; i < lights.size(); i++)
  {
    shadows = false;
    Light curLight = lights.get(i);
    PVector lightVector = new PVector(curLight.pos.x - hitInfo.hitVector.x, curLight.pos.y - hitInfo.hitVector.y, curLight.pos.z - hitInfo.hitVector.z);
    for(int k = 0; k < renderables.size(); k++)
    {
//      if(recursionLevel == 1 && k == 1)
//      {
//          println(hitInfo.hitVector);
//      }
      PVector shadowDir = lightVector;
      PVector shadowOrigin = new PVector(hitInfo.hitVector.x + shadowDir.x*0.0001, hitInfo.hitVector.y + shadowDir.y*0.0001, hitInfo.hitVector.z + shadowDir.z*0.0001);
      float lightParameter = lightVector.mag();
      Ray shadowRay = new Ray(shadowOrigin, shadowDir);
      HitInfo shadowHitInfo = renderables.get(k).intersection(shadowRay);
      if(shadowHitInfo.hit && shadowHitInfo.parameterT <= lightParameter)
      {
        shadows = true;
        k = renderables.size();
      }
    }
    
    
    
    if(!shadows)
    {
      
      //Diffuse
      
      lightVector.normalize();
      hitInfo.normal.normalize();
      
      
      
      float nDotL = max(hitInfo.normal.dot(lightVector), 0); 
    //  float nDotL = hitInfo.normal.dot(lightVector);
      PVector diffuseShading = new PVector( curRenderable.material.diffuse.x*curLight.intensity.x*nDotL, 
                                            curRenderable.material.diffuse.y*curLight.intensity.y*nDotL, 
                                            curRenderable.material.diffuse.z*curLight.intensity.z*nDotL);
                                            
                                          
      clr.x = clr.x + diffuseShading.x;
      clr.y = clr.y + diffuseShading.y;
      clr.z = clr.z + diffuseShading.z;
      
      //Specular
      PVector specularShading = computeSpecularity(curLight, renderables, curRenderable, hitInfo);
      
      
      
      
      clr.x = clr.x + specularShading.x;
      clr.y = clr.y + specularShading.y;
      clr.z = clr.z + specularShading.z;
      
      //Reflection
      
    }
//    println(curRenderable.material.reflectance);
    
  }
  
  if(curRenderable.material.reflectance > 0){
    
    PVector incidentVector = new PVector(-hitInfo.reflVector.x, -hitInfo.reflVector.y, -hitInfo.reflVector.z);
  //        println(incidentVector);
    incidentVector.normalize();
    //println(incidentVector);
    hitInfo.normal.normalize();
    //hitInfo.normal.mult(-1);
  //        println(hitInfo.normal);
    
    float dDotN = incidentVector.dot(hitInfo.normal);
    PVector secondTerm = new PVector(hitInfo.normal.x, hitInfo.normal.y, hitInfo.normal.z);
    secondTerm.mult(2*dDotN);
    PVector reflectanceVector = new PVector(incidentVector.x - secondTerm.x, incidentVector.y - secondTerm.y, incidentVector.z - secondTerm.z );
    reflectanceVector.normalize();
    
    PVector originRefl = new PVector(hitInfo.hitVector.x + 0.05*reflectanceVector.x, hitInfo.hitVector.y + 0.05*reflectanceVector.y, hitInfo.hitVector.z + 0.05*reflectanceVector.z);
    Ray reflectanceRay = new Ray(originRefl, reflectanceVector);
    int index = closestIntersection(renderables, reflectanceRay);
    
    if(index != -1 && !shadows)
    {
      //println(clr);
      Renderable reflRenderable = renderables.get(index);
      HitInfo reflHitInfo = reflRenderable.intersection(reflectanceRay);
      println(reflectanceVector);
      //reflHitInfo.normal.mult(-1);
      //println(reflRenderable.material.specular);
      PVector reflectedColor = computeColor(lights, renderables, reflRenderable, reflHitInfo, recursionLevel + 1);
      //println(reflectedColor);
      reflectedColor.mult(curRenderable.material.reflectance);
      //println(reflectedColor);
      
      //println(reflectedColor);
      clr.x = clr.x + reflectedColor.x* curRenderable.material.reflectance;
      clr.y = clr.y + reflectedColor.y* curRenderable.material.reflectance;
      clr.z = clr.z + reflectedColor.z* curRenderable.material.reflectance;
      
    }else{
      clr.x = clr.x + bgColor.x;
      clr.y = clr.y + bgColor.y;
      clr.z = clr.z + bgColor.z;
    }
  }      
  
  return clr; 
  
}

PVector computeSpecularity(Light curLight, ArrayList <Renderable> renderables, Renderable curRenderable, HitInfo hitInfo){
  PVector lightVector = new PVector(curLight.pos.x - hitInfo.hitVector.x, curLight.pos.y - hitInfo.hitVector.y, curLight.pos.z - hitInfo.hitVector.z);
  
  lightVector.normalize();
  PVector clr = new PVector(0,0,0);
  PVector halfVector = new PVector(lightVector.x, lightVector.y, lightVector.z);
  halfVector.add(hitInfo.reflVector);
  halfVector.normalize();
  float nDotH = max(hitInfo.normal.dot(halfVector), 0);
  float nDotHRaisedToP = pow(nDotH, curRenderable.material.phongExponent);
  
  PVector specularShading = new PVector(curRenderable.material.specular.x*curLight.intensity.x*nDotHRaisedToP, 
                                        curRenderable.material.specular.y*curLight.intensity.y*nDotHRaisedToP, 
                                        curRenderable.material.specular.z*curLight.intensity.z*nDotHRaisedToP);
  clr.x = specularShading.x;
  clr.y = specularShading.y;
  clr.z = specularShading.z;
  
  return clr;
}

//PVector computerReflection(){
//  
//}

class HitInfo{
  PVector reflVector;
  PVector hitVector;
  PVector normal;
  boolean hit;
  boolean triangle;
  float parameterT;
  public HitInfo( PVector rV, PVector hV,PVector n, boolean h, float pT){
    reflVector = rV;
    hitVector = hV;
    normal = n;
    hit = h;
    triangle = false;
    parameterT = pT;
  }
  public HitInfo( PVector hV,PVector n, boolean h, boolean t, float pT){
    hitVector = hV;
    normal = n;
    hit = h;
    triangle = t;
    parameterT = pT;
  }
  public HitInfo(boolean h, boolean t){
    hit = h;
    triangle = t;
  }
  
}

class Light{
  PVector pos;
  PVector intensity;
  public Light(PVector p,PVector i){
    pos = p;
    intensity = i;
  }
  
}

class Material{
  
  PVector diffuse;
  PVector ambient;
  PVector specular;
  int phongExponent;
  float reflectance;
  
  public Material(PVector d, PVector a, PVector s, int p, float refl){
    diffuse = d;
    ambient = a;
    specular = s;
    phongExponent = p;
    reflectance = refl;
  }
  public Material(Material m){
    diffuse = m.diffuse;
    ambient = m.ambient;
    specular = m.specular;
    phongExponent = m.phongExponent;
    reflectance = m.reflectance;
  }
}

class Ray{
  PVector origin;
  PVector direction;
  
  Ray(PVector o, PVector d){
    origin = o;
    direction = d;
  }
  
}

abstract class Renderable{
  
  Material material;
  
  Renderable(Material m){
    material = m;
  }

  abstract HitInfo intersection(Ray r);
  
}

class Sphere extends Renderable{
  PVector center;
  float radius;
  
  Sphere( float r, PVector c, Material m){
    super(m);
    radius = r;
    center = c;
  }
  
  HitInfo intersection(Ray r){
    HitInfo hitInfo; 
    //(d · (e − c))^2 − (d · d) ((e − c) · (e − c) − R2)
    PVector sphereVector = new PVector(r.origin.x - center.x, r.origin.y - center.y, r.origin.z - center.z);
    float firstTerm = pow(r.direction.dot(sphereVector), 2);
    
    float dirDotDir = r.direction.dot(r.direction);
    float sDotS = sphereVector.dot(sphereVector);
    float radiusSquared = pow(this.radius,2);
    float secondTerm = dirDotDir * (sDotS - radiusSquared);
    
    float discriminant = firstTerm - secondTerm;
    if(discriminant >= 0){
      float squareRootOfDiscriminant = pow(discriminant, .5);
      float negativeB = -r.direction.dot(sphereVector);
      float numerator = min(negativeB - squareRootOfDiscriminant, negativeB - squareRootOfDiscriminant);
      float parameterT = numerator/dirDotDir;
      if(parameterT > 0.0000000001)
      {
        PVector reflectionVector = new PVector(-r.direction.x, -r.direction.y, -r.direction.z );
        PVector hitPos = new PVector(r.direction.x, r.direction.y, r.direction.z);
        hitPos.mult(parameterT);
        hitPos.add(r.origin);
        
        PVector normal = new PVector((hitPos.x - center.x)/radius, (hitPos.y - center.y)/radius, (hitPos.z - center.z)/radius);  
//        PVector normal = new PVector((hitPos.x - center.x), (hitPos.y - center.y), (hitPos.z - center.z));
//        normal.mult(2);
       
//        normal.normalize();
//        reflectionVector.normalize();
        hitInfo = new HitInfo(reflectionVector, hitPos, normal, true, parameterT);
        
        return hitInfo;
      } else {
        hitInfo = new HitInfo(false, false);
        return hitInfo;
      }
    }
    else
    {
      hitInfo = new HitInfo(false, false);
      return hitInfo;
    }
  }
  
}

class Triangle extends Renderable{
  PVector pos1;
  PVector pos2;
  PVector pos3;
  int counter = 0;
  
  Triangle(PVector p1, PVector p2, PVector p3, Material m){
    super(m);
    pos1 = p1;
    pos2 = p2;
    pos3 = p3;
  }
  
  HitInfo intersection(Ray r){
   HitInfo hitInfo;
   PVector v2 = new PVector(pos2.x - pos1.x, pos2.y - pos1.y, pos2.z - pos1.z);
   PVector v1 = new PVector(pos3.x - pos1.x, pos3.y - pos1.y, pos3.z - pos1.z);
   PVector normal = v1.cross(v2);
   normal.normalize();
   this.counter++;
   
//   if(normal.dot(r.direction) == 0){
//     println(r.direction.dot(normal));
//   }
   
   if(r.direction.dot(normal) != 0)
   {
//     if(counter == (150*300 + 150)){
//       println("Mid thing");
//       println(r.direction.dot(normal));
//     }
     PVector abc = new PVector(pos1.x - pos2.x, pos1.y - pos2.y, pos1.z - pos2.z);
     PVector def = new PVector(pos1.x - pos3.x, pos1.y - pos3.y, pos1.z - pos3.z);
     PVector ghi = new PVector(r.direction.x, r.direction.y, r.direction.z);
     PVector jkl = new PVector(pos1.x - r.origin.x, pos1.y - r.origin.y, pos1.z - r.origin.z);
     
     float plane = abc.x*(def.y*ghi.z - ghi.y*def.z) + abc.y*(def.z*ghi.x - ghi.z*def.x) + abc.z*(def.x*ghi.y - ghi.x*def.y);
     float t = -(def.z*(abc.x*jkl.y - abc.y*jkl.x) + def.y*(abc.z*jkl.x - abc.x*jkl.z) + def.x*(abc.y*jkl.z - abc.z*jkl.y))/ plane;
     if(t < .0000001 )
     {
       return new HitInfo(false, true);
     }
     float gamma = (ghi.z*(abc.x*jkl.y - abc.y*jkl.x) + ghi.y*(abc.z*jkl.x - abc.x*jkl.z) + ghi.x*(abc.y*jkl.z - abc.z*jkl.y))/ plane;
     if(gamma < 0 || gamma > 1)
       return  new HitInfo(false, true);
     float beta = (jkl.x*(def.y*ghi.z - def.z * ghi.y) + jkl.y*(def.z*ghi.x - def.x * ghi.z) + jkl.z*(def.x*ghi.y - def.y * ghi.x))/ plane;
     if(beta < 0 || beta > (1 - gamma))
       return  new HitInfo(false, true);
       
     //PVector reflectionVector = new PVector(-r.direction.x, -r.direction.y, -r.direction.z );
     PVector hitPos = new PVector(r.direction.x, r.direction.y, r.direction.z);
     hitPos.mult(t);
     normal.normalize();
     //normal.mult(-1);
     PVector reflVector = new PVector(-r.direction.x, -r.direction.y, -r.direction.z);
     reflVector.normalize();
     return new HitInfo(reflVector, hitPos, normal , true, t);   
   }
   else
     return new HitInfo(false, true);
  }
  
}

void mousePressed() {
  HitInfo hitInfo = hits[mouseX][mouseY];
  if(hitInfo != null)
  {
    println("Reflection Vector" + hitInfo.reflVector);
    println("Hit Vector" + hitInfo.hitVector);
    println("X" + mouseX);
    println("Y" + mouseY);
    loadPixels();
    println(color(pixels[mouseY* 300 + mouseX]));
  }else{
    println("false");
  }
}

///////////////////////////////////////////////////////////////////////
//
// Some initializations for the scene.
//
///////////////////////////////////////////////////////////////////////
void setup() {
  size(300, 300);  
  noStroke();
  colorMode(RGB, 1.0);
  background(0, 0, 0);
  interpreter();
}

///////////////////////////////////////////////////////////////////////
//
// Draw frames.  Should leave this empty.
//
///////////////////////////////////////////////////////////////////////
void draw() {
}


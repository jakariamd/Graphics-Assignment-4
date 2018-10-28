struct Vector3
{
public:
    double x,y,z;
    Vector3() {}
    Vector3(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

class Ray
{
public:
    Vector3 start;
    Vector3 dir;
    Ray(Vector3 start, Vector3 dir)
    {
        this->start = start;
        this->dir = dir;
    }
};

class Object
{
public:
    Vector3 reference_point;
    double height, width, length;
    int shine;
    double color[3];
    double co_efficients[4];
    const double eta = 2.0;

    Object() {}
    virtual void draw() {}
    void setColor(double, double, double);
    void setShine(int shine){ this->shine = shine; };
    void setCoEfficients(double, double, double, double);
    void setRefPoint(Vector3 reference_point) { this->reference_point = reference_point; }
    void setQuadParam(double, double, double);
    virtual double intersect(Ray *r, double *current_color, int level){ return -1;}
    virtual double getIntersectingT(Ray *r){ return -1; }
};

void Object::setColor(double red, double green, double blue)
{
    color[0] = red;
    color[1] = green;
    color[2] = blue;
}

void Object::setCoEfficients(double ambient, double diffuse, double specular,double reflection)
{
    co_efficients[0] = ambient;
    co_efficients[1] = diffuse;
    co_efficients[2] = specular;
    co_efficients[3] = reflection;
}

void Object::setQuadParam(double length, double width, double height)
{
    this->height = height;
    this->width = width;
    this->length = length;
}

vector <Object*>  objects;
vector <Vector3> lights;

class Sphere:public Object
{
public:
    Sphere(Vector3 Center, double Radius)
    {
        reference_point = Center;
        length=Radius;
    }

    void draw();
    double intersect(Ray*, double *, int);
    double getIntersectingT(Ray *);
    Vector3 getNormal(Vector3);
    Vector3 getReflection(Ray*, Vector3);
    Vector3 getRefraction(Ray*, Vector3);
};

void Sphere::draw()
{
    glPushMatrix();
    {
        glTranslatef(reference_point.x,reference_point.y,reference_point.z);
        glColor3f(color[0], color[1], color[2]);
        glutSolidSphere(length,100, 100);
    }
    glPopMatrix();
}

double Sphere::getIntersectingT(Ray *r)
{
    Vector3 start;
    start.x = r->start.x - reference_point.x;
    start.y = r->start.y - reference_point.y;
    start.z = r->start.z - reference_point.z;
    double A = 1;
    double B = 2 * (r->dir.x * start.x + r->dir.y * start.y + r->dir.z * start.z);
    double C = (start.x * start.x + start.y * start.y + start.z * start.z) - length * length; 

    double D = B * B - 4 * A * C;
    if (D < 0 )
    {
        return -1;
    }

    double t1 = (-B + sqrt(D)) / float(2 * A);
    double t2 = (-B - sqrt(D)) / float(2 * A);

    if (t1 < t2) return t1;
    else return t2;
}

Vector3 Sphere::getNormal(Vector3 intersectionPoint)
{
    Vector3 normal;
    normal.x = intersectionPoint.x - reference_point.x;
    normal.y = intersectionPoint.y - reference_point.y;
    normal.z = intersectionPoint.z - reference_point.z;
    double norm = pow(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z, 0.5);
    normal.x /= norm; normal.y /= norm; normal.z /= norm;
    return normal;
}

Vector3 Sphere::getReflection(Ray* r, Vector3 intersectionPoint)
{
    Vector3 normal = getNormal(intersectionPoint);
    Vector3 reflection;
    reflection.x = r->dir.x -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.x;
    reflection.y = r->dir.y -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.y;
    reflection.z = r->dir.z -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.z;

    double norm = pow(reflection.x * reflection.x + reflection.y * reflection.y + reflection.z * reflection.z, 0.5);
    reflection.x /= norm;
    reflection.y /= norm;
    reflection.z /= norm;
    return reflection;
}

Vector3 Sphere::getRefraction(Ray* r, Vector3 normal)
{
    double N_dot_I =  normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z; 
    double k = 1.0 - eta * eta * (1.0 - N_dot_I * N_dot_I);

    Vector3 refraction(0, 0, 0);
    if (k >= 0)
    {
        refraction.x = eta * r->dir.x - (eta * N_dot_I + sqrt(k)) * normal.x;
        refraction.y = eta * r->dir.y - (eta * N_dot_I + sqrt(k)) * normal.y;
        refraction.z = eta * r->dir.z - (eta * N_dot_I + sqrt(k)) * normal.z;
        double norm = pow(refraction.x * refraction.x + refraction.y * refraction.y + refraction.z * refraction.z, 0.5);
        refraction.x /= norm; refraction.y /= norm; refraction.z /= norm;
    }

    return refraction;
}

double Sphere::intersect(Ray *r, double *current_color, int level)
{
    double t = getIntersectingT(r);
    if(t <= 0) return -1;
    if(level == 0) return t;

    Vector3 intersectionPoint;
    intersectionPoint.x = r->start.x + r->dir.x * t;
    intersectionPoint.y = r->start.y + r->dir.y * t;
    intersectionPoint.z = r->start.z + r->dir.z * t;

    for (int i=0; i<3; i++)
    {
        current_color[i] = color[i] * co_efficients[0];
    }

    Vector3 normal, reflection, refraction, start;

    for (unsigned int i = 0; i < lights.size(); ++i)
    {
        Vector3 direction ;
        direction.x = lights[i].x - intersectionPoint.x;
        direction.y = lights[i].y - intersectionPoint.y;
        direction.z = lights[i].z - intersectionPoint.z;
        double norm = pow(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z, 0.5);
        direction.x /= norm;
        direction.y /= norm;
        direction.z /= norm;

        start.x = intersectionPoint.x + direction.x ;
        start.y = intersectionPoint.y + direction.y ;
        start.z = intersectionPoint.z + direction.z ;

        Ray L(start, direction);
        bool is_obscured = false;
        for (unsigned int j = 0; j < objects.size(); ++j)
        {
            double t_obscure = objects[j]->getIntersectingT(&L);
            if (t_obscure >= 0 && t_obscure < norm)
            {
                is_obscured = true;
                break;
            }
        }
        if (!is_obscured)
        {
            normal = getNormal(intersectionPoint);
            double lambert_value = (L.dir.x*normal.x+L.dir.y*normal.y+L.dir.z*normal.z) * co_efficients[1] ;
            if(lambert_value < 0) lambert_value = 0;

            reflection = getReflection(r, intersectionPoint);
            double phong_value = pow(reflection.x * r->dir.x + reflection.y * r->dir.y + reflection.z*r->dir.z, shine) * co_efficients[2];
            if (phong_value < 0) phong_value = 0;

            for (int k=0; k<3; k++)
            {
                current_color[k] += (lambert_value + phong_value) * color[k];
            }
        }

        if(level < level_of_recursion)
        {
            start.x = intersectionPoint.x + reflection.x ;
            start.y = intersectionPoint.y + reflection.y ;
            start.z = intersectionPoint.z + reflection.z ;

            Ray reflectionRay(start, reflection);
            double nearest = -1;
            double dummyColorAt[3];
            double t_min = INT_MAX;
            for (unsigned int k = 0; k < objects.size(); ++k)
            {
                double t_near = objects[k]->getIntersectingT(&reflectionRay);
                if(t_near <= 0) continue;
                if (t_near < t_min)
                {
                    t_min = t_near;
                    nearest = k;
                }
            }
            if (nearest != -1)
            {
                double t = objects[nearest]->intersect(&reflectionRay, dummyColorAt, level + 1);
                for (int k=0; k<3; k++)
                {
                    current_color[k] += dummyColorAt[k] * co_efficients[3];
                }
            }


            refraction = getRefraction(r, normal);
            start.x = intersectionPoint.x + refraction.x ;
            start.y = intersectionPoint.y + refraction.y ;
            start.z = intersectionPoint.z + refraction.z ;

            Ray refractionRay(start, refraction);
            nearest = -1;
            dummyColorAt[0] = 0; dummyColorAt[1] = 0 ; dummyColorAt[2] =  0;
            t_min = INT_MAX;
            for (unsigned int k = 0; k < objects.size(); ++k)
            {
                double t_near = objects[k]->getIntersectingT(&refractionRay);
                if(t_near <= 0) continue;
                if (t_near < t_min)
                {
                    t_min = t_near;
                    nearest = k;
                }
            }
            if (nearest != -1)
            {
                double t = objects[nearest]->intersect(&refractionRay, dummyColorAt, level + 1);
                for (int k=0; k<3; k++)
                {
                    current_color[k] += dummyColorAt[k] * eta;
                }
            }
        }


        for (int i = 0; i < 3; ++i)
        {
            if (current_color[i] < 0) current_color[i] = 0;
            if (current_color[i] > 1) current_color[i] = 1;
        }

    }

    return t;
}

class Floor: public Object
{
public:
    bitmap_image texImage;
    Floor(double FloorWidth, double TileWidth)
    {
        reference_point = Vector3(-FloorWidth/2, -FloorWidth/2,0);
        length=TileWidth;
        texImage = bitmap_image("flower.bmp");

    }

    void draw();
    double getIntersectingT(Ray *);
    Vector3 getNormal();
    double intersect(Ray*, double *, int);
    Vector3 getReflection(Ray*);
};

void Floor::draw()
{
    int tile_count = -reference_point.x * 2 / length;
    for (int i = 0; i < tile_count; ++i)
    {
        for (int j = 0; j < tile_count; ++j)
        {
            int tile_color = (i+j) % 2;
            glColor3f(tile_color, tile_color, tile_color);
            glBegin(GL_QUADS);
            {
                glVertex3f(reference_point.x + i * length, reference_point.y + j * length, 0);
                glVertex3f(reference_point.x + (i+1) * length, reference_point.y + j * length, 0);
                glVertex3f(reference_point.x + (i+1) * length, reference_point.y + (j+1) * length, 0);
                glVertex3f(reference_point.x + i * length, reference_point.y + (j+1) * length, 0);
            }
            glEnd();
        }
    }
}

Vector3 Floor::getNormal()
{
    Vector3 normal(0, 0, 1);
    return normal;
}

double Floor::getIntersectingT(Ray* r)
{
    return -r->start.z / r->dir.z;
}

Vector3 Floor::getReflection(Ray* r)
{
    Vector3 reflection(r->dir.x, r->dir.y, -r->dir.z);
    return reflection;
}

double Floor::intersect(Ray *r, double *current_color, int level)
{
    double t = getIntersectingT(r);

    Vector3 intersectionPoint;
    intersectionPoint.x = r->start.x + r->dir.x * t;
    intersectionPoint.y = r->start.y + r->dir.y * t;
    intersectionPoint.z = r->start.z + r->dir.z * t;

    if (intersectionPoint.x < reference_point.x || intersectionPoint.y < reference_point.y ||
            intersectionPoint.x > (-reference_point.x) || intersectionPoint.y > (-reference_point.y))
    {
        return -1;
    }
    if(t <= 0) return -1;
    if(level == 0) return t;


    unsigned char red, green, blue;
    int x = (intersectionPoint.x + abs(reference_point.x) - int((intersectionPoint.x - reference_point.x)/length) * length) *  texImage.width()/length;
    int y = (intersectionPoint.y + abs(reference_point.y) - int((intersectionPoint.y - reference_point.y)/length) * length) *  texImage.height()/length;
    texImage.get_pixel(x, y, red, green, blue);
    double texColor[] = {red, green, blue};

    int col = (int((intersectionPoint.x - reference_point.x)/length) + int((intersectionPoint.y - reference_point.y)/length)) % 2;
    for (int i=0; i<3; i++)
    {
        color[i] = col * texColor[i] / 255.0;
        current_color[i] = color[i] * co_efficients[0] ;
        //color[i] = col;
        //current_color[i] = col * co_efficients[0] ;
    }

    Vector3 normal(0, 0, 1), start, reflection = getReflection(r);
    for (unsigned int i = 0; i < lights.size(); ++i)
    {
        Vector3 direction ;
        direction.x = lights[i].x - intersectionPoint.x;
        direction.y = lights[i].y - intersectionPoint.y;
        direction.z = lights[i].z - intersectionPoint.z;
        double norm = pow(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z, 0.5);
        direction.x /= norm;
        direction.y /= norm;
        direction.z /= norm;

        start.x = intersectionPoint.x + direction.x ;
        start.y = intersectionPoint.y + direction.y ;
        start.z = intersectionPoint.z + direction.z ;

        Ray L(start, direction);


        bool is_obscured = false;
        for (unsigned int j = 0; j < objects.size(); ++j)
        {
            double t_obscure = objects[j]->getIntersectingT(&L);
            if (t_obscure >= 0 && t_obscure < norm)
            {
                is_obscured = true;
                break;
            }
        }

        if (!is_obscured)
        {
            double lambert_value = (L.dir.x*normal.x + L.dir.y*normal.y + L.dir.z*normal.z) * co_efficients[1] ;
            if(lambert_value < 0) lambert_value = 0;

            double phong_value = pow(reflection.x * r->dir.x + reflection.y * r->dir.y + reflection.z*r->dir.z, shine) * co_efficients[2];
            if (phong_value < 0) phong_value = 0;

            for (int k=0; k<3; k++)
            {
                current_color[k] += 1.0 * (lambert_value + phong_value) * color[k];
            }
        }


        if(level < level_of_recursion)
        {
            start.x = intersectionPoint.x + reflection.x ;
            start.y = intersectionPoint.y + reflection.y ;
            start.z = intersectionPoint.z + reflection.z ;

            Ray reflectionRay(start, reflection);
            double nearest = -1;
            double dummyColorAt[3];
            double t_min = INT_MAX;
            for (unsigned int k = 0; k < objects.size(); ++k)
            {
                double t_near = objects[k]->getIntersectingT(&reflectionRay);
                if(t_near < 0) continue;
                if (t_near < t_min)
                {
                    t_min = t_near;
                    nearest = k;
                }
            }
            if (nearest != -1)
            {
                double t = objects[nearest]->intersect(&reflectionRay, dummyColorAt, level + 1);
                for (int k=0; k<3; k++)
                {
                    current_color[k] += dummyColorAt[k] * co_efficients[3];
                }
            }
        }

        for (int i = 0; i < 3; ++i)
        {
            if (current_color[i] < 0) current_color[i] = 0;
            if (current_color[i] > 1) current_color[i] = 1;
        }

    }

    return t;
}



class Triangle:public Object
{
    Vector3 point1, point2, point3;
public:
    Triangle(Vector3 point1, Vector3 point2, Vector3 point3)
    {
        this->point1 = point1;
        this->point2 = point2;
        this->point3 = point3;
    }

    void draw();
    double getIntersectingT(Ray *);
    Vector3 getNormal();
    double intersect(Ray*, double *, int);
    Vector3 getReflection(Ray*, Vector3 normal);

};

void Triangle::draw()
{
    glColor3f(color[0],color[1],color[2]);
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(point1.x, point1.y, point1.z);
        glVertex3f(point2.x, point2.y, point2.z);
        glVertex3f(point3.x, point3.y, point3.z);
    }
    glEnd();
}

Vector3 Triangle::getNormal()
{
    Vector3 u, v, normal;
    u.x = point2.x - point1.x;
    u.y = point2.y - point1.y;
    u.z = point2.z - point1.z;
    v.x = point3.x - point1.x;
    v.y = point3.y - point1.y;
    v.z = point3.z - point1.z;

    normal.x = u.y*v.z - u.z*v.y;
    normal.y = u.z*v.x - u.x*v.z;
    normal.z = u.x*v.y - u.y*v.x;
    double norm = pow(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z, 0.5);
    normal.x /= norm;
    normal.y /= norm;
    normal.z /= norm;
    return normal;
}


double Triangle::getIntersectingT(Ray* r)
{
    const float EPSILON = 0.0000001;
    Vector3 edge1, edge2, h, s, q;
    float a,f,u,v;
    edge1.x = point2.x - point1.x;
    edge1.y = point2.y - point1.y;
    edge1.z = point2.z - point1.z;
    edge2.x = point3.x - point1.x;
    edge2.y = point3.y - point1.y;
    edge2.z = point3.z - point1.z;

    h.x = r->dir.y * edge2.z - r->dir.z * edge2.y;
    h.y = r->dir.z * edge2.x - r->dir.x * edge2.z;
    h.z = r->dir.x * edge2.y - r->dir.y * edge2.x;

    a = h.x * edge1.x + h.y * edge1.y + h.z * edge1.z;
    if (a > -EPSILON && a < EPSILON) return -1;

    f = 1/a;
    s.x = r->start.x - point1.x;
    s.y = r->start.y - point1.y;
    s.z = r->start.z - point1.z;
    u = f * (s.x * h.x + s.y * h.y + s.z * h.z);
    if (u < 0.0 || u > 1.0) return -1;

    q.x = s.y * edge1.z - s.z*edge1.y;
    q.y = s.z * edge1.x - s.x*edge1.z;
    q.z = s.x * edge1.y - s.y*edge1.x;

    v = f * (r->dir.x * q.x + r->dir.y * q.y + r->dir.z * q.z);
    if (v < 0.0 || u + v > 1.0) return -1;

    double t = f * (edge2.x * q.x + edge2.y * q.y + edge2.z * q.z);
    if (t > EPSILON)
    {
        return t;
    }
    else return -1;
}


Vector3 Triangle::getReflection(Ray* r, Vector3 normal)
{
    Vector3 reflection;
    reflection.x = r->dir.x -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.x;
    reflection.y = r->dir.y -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.y;
    reflection.z = r->dir.z -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.z;

    double norm = pow(reflection.x * reflection.x + reflection.y * reflection.y + reflection.z * reflection.z, 0.5);
    reflection.x /= norm;
    reflection.y /= norm;
    reflection.z /= norm;
    return reflection;
}

double Triangle::intersect(Ray *r, double *current_color, int level)
{
    double t = getIntersectingT(r);
    if(t <= 0) return -1;
    if(level == 0) return t;

    Vector3 intersectionPoint;
    intersectionPoint.x = r->start.x + r->dir.x * t;
    intersectionPoint.y = r->start.y + r->dir.y * t;
    intersectionPoint.z = r->start.z + r->dir.z * t;

    for (int i=0; i<3; i++)
    {
        current_color[i] = color[i] * co_efficients[0];
    }

    Vector3 normal, start;
    for (unsigned int i = 0; i < lights.size(); ++i)
    {
        Vector3 direction ;
        direction.x = lights[i].x - intersectionPoint.x; 
        direction.y = lights[i].y - intersectionPoint.y;
        direction.z = lights[i].z - intersectionPoint.z;
        double norm = pow(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z, 0.5);
        direction.x /= norm;
        direction.y /= norm;
        direction.z /= norm;

        start.x = intersectionPoint.x + direction.x ;
        start.y = intersectionPoint.y + direction.y ;
        start.z = intersectionPoint.z + direction.z ;

        normal = getNormal();
        if ((direction.x * normal.x + direction.y * normal.y  + direction.z * normal.z) > 0)
        {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }

        Vector3 reflection = getReflection(r, normal);

        Ray L(start, direction);
        bool is_obscured = false;
        for (unsigned int j = 0; j < objects.size(); ++j)
        {
            double t_obscure = objects[j]->getIntersectingT(&L);
            if (t_obscure >= 0 && t_obscure < norm)
            {
                is_obscured = true;
                break;
            }
        }
        if (!is_obscured)
        {
            double lambert_value = (L.dir.x*normal.x + L.dir.y*normal.y + L.dir.z*normal.z) * co_efficients[1] ;
            if(lambert_value < 0) lambert_value = 0;

            double phong_value = pow((reflection.x * r->dir.x + reflection.y * r->dir.y + reflection.z*r->dir.z),
                                     shine) * co_efficients[2];
            if (phong_value < 0) phong_value = 0;

            for (int k=0; k<3; k++)
            {
                current_color[k] += (lambert_value + phong_value) * color[k];
            }
        }


        if(level < level_of_recursion)
        {
            start.x = intersectionPoint.x + reflection.x ;
            start.y = intersectionPoint.y + reflection.y ;
            start.z = intersectionPoint.z + reflection.z ;

            Ray reflectionRay(start, reflection);
            double nearest = -1;
            double dummyColorAt[3];
            double t_min = INT_MAX;
            for (unsigned int k = 0; k < objects.size(); ++k)
            {
                double t_near = objects[k]->getIntersectingT(&reflectionRay);
                if(t_near < 0) continue;
                if (t_near < t_min)
                {
                    t_min = t_near;
                    nearest = k;
                }
            }
            if (nearest != -1)
            {
                double t = objects[nearest]->intersect(&reflectionRay, dummyColorAt, level + 1);
                if (t >= 0)
                {
                    for (int k=0; k<3; k++)
                    {
                        current_color[k] += dummyColorAt[k] * co_efficients[3];
                    }
                }
            }
        }

        for (int i = 0; i < 3; ++i)
        {
            if (current_color[i] < 0) current_color[i] = 0;
            if (current_color[i] > 1) current_color[i] = 1;
        }
    }

    return t;
}




class General:public Object
{
    double A, B, C, D, E, F, G, H, I, J;
public:
    General(double A, double B, double C, double D, double E, double F, double G, double H, double I, double J)
    {
        this->A = A; this->B = B; this->C = C; this->D = D; this->E = E;
        this->F = F; this->G = G; this->H = H; this->I = I; this->J = J;
    }

    void draw() {}
    Vector3 getNormal(Vector3);
    double getIntersectingT(Ray *);
    double intersect(Ray*, double *, int);
    Vector3 getReflection(Ray*, Vector3);
    bool checkLimit(Vector3 *);
};

Vector3 General::getNormal(Vector3 intersectionPoint){
    Vector3 normal;
    normal.x = 2 * A * intersectionPoint.x + D * intersectionPoint.y + F * intersectionPoint.z  + G;
    normal.y = 2 * B * intersectionPoint.y + D * intersectionPoint.x + E * intersectionPoint.z  + H;
    normal.z = 2 * C * intersectionPoint.z + E * intersectionPoint.y + F * intersectionPoint.x  + I;
    double norm = pow(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z, 0.5);
    normal.x /= norm;
    normal.y /= norm;
    normal.z /= norm;
    return normal;
}

bool General::checkLimit(Vector3* intPiont){
    bool c1 = length > 0 && ( reference_point.x > intPiont->x || intPiont->x > (reference_point.x + length));
    bool c2 = width > 0 && (reference_point.y > intPiont->y  || intPiont->y > (reference_point.y + width));  
    bool c3 = height > 0 && (reference_point.z > intPiont->z || intPiont->z > (reference_point.z + height));
    return (c1||c2||c3); 
}
double General::getIntersectingT(Ray* r){
    double a = A * pow(r->dir.x,2) + B * pow(r->dir.y,2) + C * pow(r->dir.z,2) +
               D * r->dir.x * r->dir.y + E * r->dir.y * r->dir.z + F * r->dir.z * r->dir.x;

    double b = 2 * (A * r->start.x * r->dir.x + B * r->start.y * r->dir.y + C * r->start.z * r->dir.z) +
               D * (r->start.x * r->dir.y + r->dir.x * r->start.y) + E * (r->start.y * r->dir.z + r->dir.y * r->start.z) +
               F * (r->start.z * r->dir.x + r->dir.z * r->start.x) + G * r->dir.x + H * r->dir.y + I * r->dir.z;

    double c = A * pow(r->start.x,2) + B * pow(r->start.y,2) + C * pow(r->start.z,2) +
               D * r->start.x * r->start.y + E * r->start.y * r->start.z + F * r->start.z * r->start.x +
               G * r->start.x + H * r->start.y + I * r->start.z + J;

    double d = b*b - 4*a*c;
    if (d < 0 )
    {
        return -1;
    }

    double t1 = (-b + sqrt(d)) / float(2 * a);
    double t2 = (-b - sqrt(d)) / float(2 * a);

    Vector3 intPiont1, intPiont2;
    intPiont1.x = r->start.x + r->dir.x * t1; intPiont1.y = r->start.y + r->dir.y * t1; intPiont1.z = r->start.z + r->dir.z * t1;
    intPiont2.x = r->start.x + r->dir.x * t2; intPiont2.y = r->start.y + r->dir.y * t2; intPiont2.z = r->start.z + r->dir.z * t2;

    bool cond1 = checkLimit(&intPiont1);
    bool cond2 = checkLimit(&intPiont2);

    if (cond1 && cond2) return -1;
    if (!cond1) return t1;
    if (!cond2) return t2;
    if (t1 < t2) return t1;
    else return t2;


}

Vector3 General::getReflection(Ray* r, Vector3 normal)
{
    Vector3 reflection;
    reflection.x = r->dir.x -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.x;
    reflection.y = r->dir.y -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.y;
    reflection.z = r->dir.z -  2 * (normal.x * r->dir.x + normal.y * r->dir.y + normal.z * r->dir.z) * normal.z;

    double norm = pow(reflection.x * reflection.x + reflection.y * reflection.y + reflection.z * reflection.z, 0.5);
    reflection.x /= norm;
    reflection.y /= norm;
    reflection.z /= norm;
    return reflection;
}


double General::intersect(Ray *r, double *current_color, int level)
{
    double t = getIntersectingT(r);
    if(t <= 0) return -1;
    if(level == 0) return t;

    Vector3 intersectionPoint;
    intersectionPoint.x = r->start.x + r->dir.x * t;
    intersectionPoint.y = r->start.y + r->dir.y * t;
    intersectionPoint.z = r->start.z + r->dir.z * t;

    for (int i=0; i<3; i++)
    {
        current_color[i] = color[i] * co_efficients[0];
    }

    for (unsigned int i = 0; i < lights.size(); ++i)
    {
        Vector3 direction ;
        direction.x = lights[i].x - intersectionPoint.x; 
        direction.y = lights[i].y - intersectionPoint.y;
        direction.z = lights[i].z - intersectionPoint.z;
        double norm = pow(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z, 0.5);
        direction.x /= norm;
        direction.y /= norm;
        direction.z /= norm;

        Vector3 start;
        start.x = intersectionPoint.x + direction.x ;
        start.y = intersectionPoint.y + direction.y ;
        start.z = intersectionPoint.z + direction.z ;

        Vector3 normal = getNormal(intersectionPoint);
        Vector3 reflection = getReflection(r, normal);

        Ray L(start, direction);
        bool is_obscured = false;

        for (unsigned int j = 0; j < objects.size(); ++j)
        {
            double t_obscure = objects[j]->getIntersectingT(&L);
            if (t_obscure >= 0 && t_obscure < norm)
            {
                is_obscured = true;
                break;
            }
        }


        if (!is_obscured)
        {
            double lambert_value = (L.dir.x*normal.x+L.dir.y*normal.y+L.dir.z*normal.z) * co_efficients[1] ;
            if(lambert_value < 0) lambert_value = 0;

            double phong_value = pow((reflection.x * r->dir.x + reflection.y * r->dir.y + reflection.z*r->dir.z),
                                     shine) * co_efficients[2];
            if (phong_value < 0) phong_value = 0;

            for (int k=0; k<3; k++)
            {
                current_color[k] += (lambert_value + phong_value) * color[k];
            }
        }

        if(level < level_of_recursion)
        {
            start.x = intersectionPoint.x + reflection.x ;
            start.y = intersectionPoint.y + reflection.y ;
            start.z = intersectionPoint.z + reflection.z ;

            Ray reflectionRay(start, reflection);
            double nearest = -1;
            double dummyColorAt[3];
            double t_min = INT_MAX;
            for (unsigned int k = 0; k < objects.size(); ++k)
            {
                double t_near = objects[k]->getIntersectingT(&reflectionRay);
                if(t_near < 0) continue;
                if (t_near < t_min)
                {
                    t_min = t_near;
                    nearest = k;
                }
            }
            if (nearest != -1)
            {
                double t = objects[nearest]->intersect(&reflectionRay, dummyColorAt, level + 1);
                for (int k=0; k<3; k++)
                {
                    current_color[k] += dummyColorAt[k] * co_efficients[3];
                }
            }
        }

        for (int i = 0; i < 3; ++i)
        {
            if (current_color[i] < 0) current_color[i] = 0;
            if (current_color[i] > 1) current_color[i] = 1;
        }

    }

    return t;
}
